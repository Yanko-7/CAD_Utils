use plotters::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open("data.txt")?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let points = read_points(&mut lines)?;
    let boundary_indices = read_boundary_indices(&mut lines)?;
    let triangulation_edges = read_triangulation_edges(&mut lines)?;
    let constraint_edges = read_constraint_edges(&mut lines)?;

    draw_plot(&points, &boundary_indices, &triangulation_edges, &constraint_edges)?;

    Ok(())
}

fn read_points(lines: &mut std::io::Lines<BufReader<File>>) -> Result<Vec<(f64, f64)>, Box<dyn std::error::Error>> {
    let mut points = Vec::new();
    while let Some(line) = lines.next() {
        let line = line?;
        if line == "Boundary Indices:" {
            break;
        }
        if line == "Points:" {
            continue;
        }
        let coords: Vec<&str> = line.split_whitespace().collect();
        let x: f64 = coords[0].parse()?;
        let y: f64 = coords[1].parse()?;
        points.push((x, y));
    }
    Ok(points)
}

fn read_boundary_indices(lines: &mut std::io::Lines<BufReader<File>>) -> Result<Vec<usize>, Box<dyn std::error::Error>> {
    let mut boundary_indices = Vec::new();
    while let Some(line) = lines.next() {
        let line = line?;
        if line == "Triangulation Edges:" {
            break;
        }
        let idx: usize = line.parse()?;
        boundary_indices.push(idx);
    }
    Ok(boundary_indices)
}

fn read_triangulation_edges(lines: &mut std::io::Lines<BufReader<File>>) -> Result<Vec<(usize, usize)>, Box<dyn std::error::Error>> {
    let mut edges = Vec::new();
    while let Some(line) = lines.next() {
        let line = line?;
        if line == "Constraint Edges:" {
            break;
        }
        let indices: Vec<&str> = line.split_whitespace().collect();
        let start: usize = indices[0].parse()?;
        let end: usize = indices[1].parse()?;
        edges.push((start, end));
    }
    Ok(edges)
}

fn read_constraint_edges(lines: &mut std::io::Lines<BufReader<File>>) -> Result<Vec<(usize, usize)>, Box<dyn std::error::Error>> {
    let mut edges = Vec::new();
    while let Some(line) = lines.next() {
        let line = line?;
        let indices: Vec<&str> = line.split_whitespace().collect();
        let start: usize = indices[0].parse()?;
        let end: usize = indices[1].parse()?;
        edges.push((start, end));
    }
    Ok(edges)
}

fn draw_plot(
    points: &[(f64, f64)],
    boundary_indices: &[usize],
    triangulation_edges: &[(usize, usize)],
    constraint_edges: &[(usize, usize)],
) -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new("plot.svg", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(10, 10, 10, 10);
    let mut chart = ChartBuilder::on(&root)
        .caption("Points and Edges", ("sans-serif", 30).into_font())
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(-1.0..2.0, -1.0..2.0)?;
    chart.configure_mesh().draw()?;

    chart.draw_series(points.iter().map(|(x, y)| Circle::new((*x, *y), 1, BLACK.filled())))?;
    chart.draw_series(boundary_indices.iter().map(|&idx| {
        let (x, y) = points[idx];
        Circle::new((x, y), 1, RED.filled())
    }))?;
    chart.draw_series(triangulation_edges.iter().map(|&(start, end)| {
        let (x0, y0) = points[start];
        let (x1, y1) = points[end];
        PathElement::new(vec![(x0, y0), (x1, y1)], &BLACK)
    }))?;
    chart.draw_series(constraint_edges.iter().map(|&(start, end)| {
        let (x0, y0) = points[start];
        let (x1, y1) = points[end];
        PathElement::new(vec![(x0, y0), (x1, y1)], &BLUE)
    }))?;

    root.present()?;
    Ok(())
}
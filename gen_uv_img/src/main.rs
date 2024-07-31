use plotters::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open("C:/Users/yanko/GME/build/src_acis/test/data.txt")?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let u = read_interval(&mut lines, "Uinterval:")?;
    let v = read_interval(&mut lines, "Vinterval:")?;
    let points = read_points(&mut lines)?;
    let boundary_indices = read_boundary_indices(&mut lines)?;
    let sample_indices = read_sample_indices(&mut lines)?;
    let triangulation_edges = read_triangulation_edges(&mut lines)?;
    let constraint_edges = read_constraint_edges(&mut lines)?;

    draw_plot(&points, &boundary_indices, &sample_indices, &triangulation_edges, &constraint_edges, u, v)?;

    Ok(())
}
fn read_interval(lines: &mut std::io::Lines<BufReader<File>>, label: &str) -> Result<(f64, f64), Box<dyn std::error::Error>> {
    while let Some(line) = lines.next() {
        let line = line?;
        if line == label {
            let interval_line = lines.next().ok_or("Missing interval line")??;
            let coords: Vec<&str> = interval_line.split_whitespace().collect();
            let start: f64 = coords[0].parse()?;
            let end: f64 = coords[1].parse()?;
            return Ok((start, end));
        }
    }
    Err("Interval not found".into())
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
        if line == "Sample Indices:" {
            break;
        }
        let idx: usize = line.parse()?;
        boundary_indices.push(idx);
    }
    Ok(boundary_indices)
}
fn read_sample_indices(lines: &mut std::io::Lines<BufReader<File>>) -> Result<Vec<usize>, Box<dyn std::error::Error>> {
    let mut sample_indices = Vec::new();
    while let Some(line) = lines.next() {
        let line = line?;
        if line == "Triangulation Edges:" {
            break;
        }
        let idx: usize = line.parse()?;
        sample_indices.push(idx);
    }
    Ok(sample_indices)
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
    sample_indices: &[usize],
    triangulation_edges: &[(usize, usize)],
    constraint_edges: &[(usize, usize)],
    u: (f64, f64),
    v: (f64, f64),
    //参数空间
) -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new("plot.svg", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(10, 10, 10, 10);
    let mut chart = ChartBuilder::on(&root)
        .caption("Parameter Space", ("sans-serif", 20).into_font())
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(u.0..u.1, v.0..v.1)?;
    chart.configure_mesh().draw()?;

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
    chart.draw_series(points.iter().map(|(x, y)| Circle::new((*x, *y), 1, GREEN.filled())))?.label("Intersection")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart.draw_series(boundary_indices.iter().map(|&idx| {
        let (x, y) = points[idx];
        Circle::new((x, y), 1, YELLOW.filled())
    }))?.label("boundary point")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &YELLOW));
    
    chart.draw_series(sample_indices.iter().map(|&idx| {
        let (x, y) = points[idx];
        Circle::new((x, y), 1, RED.filled())
    }))?.label("sampling point")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    // 配置标签样式
    // 采样点
    chart.configure_series_labels()
    .background_style(&WHITE.mix(0.8))
    .border_style(&BLACK)
    .draw()?;
 
    root.present()?;
    Ok(())
}
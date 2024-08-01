import os

# Constants
FOLDER_PATH = 'acis2gme_test'
DISTANCE_THRESHOLD = 4


def process_line(content: str):
    """Process a line to check for specific patterns and modify it."""
    if "#include" in content:
        return content, False
    modified_content = ''
    flag = False
    i = 0
    while i < len(content):
        if content[i:i + 4] == 'gme_':
            if content[i + 4:i + 9] != 'facet':
                flag = True
                i += 4
                continue
        modified_content += content[i]
        i += 1
    return modified_content, flag


def process_file(file_path: str):
    """Process a single file and apply necessary modifications."""
    new_lines = []
    flag = False
    begin_idx = -DISTANCE_THRESHOLD
    end_idx = -DISTANCE_THRESHOLD

    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

        for i, line in enumerate(lines):
            modified_line, is_gme = process_line(line)
            if is_gme:
                if not flag:
                    flag = True
                    begin_idx = i
                end_idx = i
            else:
                if i - end_idx >= DISTANCE_THRESHOLD:
                    if flag:
                        flag = False
                        new_lines.append("#ifdef FACETER_USE_ACIS\n")
                        for j in range(begin_idx, end_idx + 1):
                            new_str, _ = process_line(lines[j])
                            new_lines.append(new_str)
                        new_lines.append("#else\n")
                        for j in range(begin_idx, end_idx + 1):
                            new_lines.append(lines[j])
                        new_lines.append("#endif\n")
                        for j in range(end_idx + 1, i + 1):
                            new_lines.append(lines[j])
                    else:
                        new_lines.append(line)

        if flag:
            new_lines.append("#ifdef FACETER_USE_ACIS\n")
            for j in range(begin_idx, end_idx + 1):
                new_str, _ = process_line(lines[j])
                new_lines.append(new_str)
            new_lines.append("#else\n")
            for j in range(begin_idx, end_idx + 1):
                new_lines.append(lines[j])
            new_lines.append("#endif\n")
            for j in range(end_idx + 1, i + 1):
                new_lines.append(lines[j])

    return new_lines


def process_all_files_in_folder(folder_path: str):
    """Process all .cpp files in the given folder."""
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)

        if os.path.isfile(file_path) and filename.endswith('.cpp'):
            modified_content = process_file(file_path)
            with open(file_path, 'w', encoding='utf-8') as file:
                file.writelines(modified_content)


if __name__ == '__main__':
    process_all_files_in_folder(FOLDER_PATH)

import os
import re

folder_path = 'gme2acis_test'
output_folder_path = 'gme2acis_test_output'
DISTANCE = 4


def remove_substrings(input_string):
    """Remove occurrences of 'gme' substrings from the input string."""
    res = input_string.replace("\"gme\",", "")
    res = res.replace("AcisOptions", "GmeOptions")
    res = res.replace("ACIS_DELETE ", "GME_DELETE ")
    res = res.replace("ACIS_NEW ", "GME_NEW ")
    res = res.replace("ACIS_DELETE", "GME_DELETE")
    res = res.replace("API_BEGIN", "__GME_API_BEGIN")
    res = res.replace("API_END", "__GME_API_END")
    return res

def process_line(content):
    """
    Process a line of content.
    - Removes 'gme_' prefixes unless followed by 'facet'.
    - Removes 'gme' substrings.
    """
    if '#include' in content:
        return content, False

    flag = False
    processed_content = ''
    content_length = len(content)
    index = 0

    while index < content_length:
        if index + 4 < content_length and content[index:index + 4] == 'gme_':
            if index + 9 <= content_length and content[index + 4: index + 9] != 'facet':
                flag = True
                index += 4
                continue
        processed_content += content[index]
        index += 1

    return remove_substrings(processed_content), flag


def process_file(file_path):
    """
    Process a single file:
    - Removes lines between '#ifdef FACETER_USE' and '#else'.
    - Keeps lines between '#else' and '#endif'.
    - Calls process_line function to process each line.
    """
    processed_lines = []
    status = 0

    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    for line in lines:
        new_line, _ = process_line(line)
        processed_lines.append(new_line)

    return processed_lines


def process_files_in_directory(directory_path, output_directory):
    """
    Process all .cpp files in the given directory:
    - Calls process_file function for each file.
    - Saves processed content to a new file, removing 'gme_' prefix from filename.
    """
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)

        if os.path.isfile(file_path) and filename.endswith('.cpp'):
            processed_content = process_file(file_path)
            new_filename = filename
            if filename.startswith("gme_"):
                new_filename = filename[len("gme_"):]
            new_file_path = os.path.join(output_directory, new_filename)
            with open(new_file_path, 'w', encoding='utf-8') as file:
                file.writelines(processed_content)


if __name__ == '__main__':
    process_files_in_directory(folder_path,output_folder_path )

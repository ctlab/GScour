import re
"""
script for collect files in broken_list to relaunch paml or for other processing
a temporary solution as long as the list is not in 
shared variables between processes
"""
broken_list = list()
with open('/path/to/paml_one_ratio.log', 'r') as f:
    for line in f:
        if "Null size result file number" in line or "Not null return code" in line:
            file_number = (re.search(r"number\s(\d+)", line)).group(1)
            broken_list.append(file_number)

print(broken_list)
print(len(broken_list))
import re
import os
import sys

# Usage: python fortran_deps.py src_files.txt build_dir
# src_files.txt contains list of Fortran source files, one per line
# build_dir is the directory where .o files are placed (to rewrite .o paths)

def parse_modules(filename):
    """Parse Fortran file for modules defined and modules used."""
    defined = set()
    used = set()
    module_re = re.compile(r'^\s*module\s+(\w+)', re.I)
    use_re = re.compile(r'^\s*use\s+(\w+)', re.I)
    ignore_re = re.compile(r'^\s*module procedure', re.I)

    with open(filename, 'r') as f:
        for line in f:
            # Skip "module procedure" lines (not module definitions)
            if ignore_re.match(line):
                continue

            m = module_re.match(line)
            if m:
                modname = m.group(1).lower()
                # Ignore 'module' statements inside program or subroutine if any
                if modname != 'procedure':
                    defined.add(modname)
                continue

            u = use_re.match(line)
            if u:
                used.add(u.group(1).lower())

    return defined, used

def main():
    if len(sys.argv) != 3:
        print("Usage: python fortran_deps.py src_files.txt build_dir")
        sys.exit(1)

    src_list_file = sys.argv[1]
    build_dir = sys.argv[2].rstrip('/')

    # Read source files
    with open(src_list_file) as f:
        src_files = [line.strip() for line in f if line.strip()]

    # Map module name -> source file defining it
    module_to_file = {}

    # Map source file -> modules it defines and uses
    file_modules = {}

    for src in src_files:
        defined, used = parse_modules(src)
        file_modules[src] = {'defined': defined, 'used': used}
        for m in defined:
            if m in module_to_file:
                print(f"Warning: module {m} defined in multiple files: {module_to_file[m]} and {src}")
            module_to_file[m] = src

    # Build dependency map: for each source file, list files it depends on (via used modules)
    file_deps = {}

    for src, mods in file_modules.items():
        deps = set()
        for used_mod in mods['used']:
            if used_mod in module_to_file and module_to_file[used_mod] != src:
                deps.add(module_to_file[used_mod])
        file_deps[src] = deps

    # Now print Makefile style dependencies:
    # Convert src paths to object paths in build_dir
    def obj_file(src):
        # Replace src/ prefix with build_dir/
        # handle subdirs properly
        if src.startswith('src/'):
            rel = src[4:]
        else:
            rel = src
        return os.path.join(build_dir, os.path.splitext(rel)[0] + '.o')

    for src in src_files:
        obj = obj_file(src)
        dep_objs = [obj_file(dep) for dep in sorted(file_deps[src])]
        deps_str = ' '.join(dep_objs)
        print(f"{obj}: {deps_str}")

if __name__ == "__main__":
    main()

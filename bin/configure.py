import os
import sys
import shutil

def parse_arguments():
    """
    Parse command-line arguments in the form +key=value.
    Returns a dictionary of arguments with default values.

    Defaults:
        simulation: 'sedov'
        hydro: 'rk2'
        make: 'linux'
        objdir: 'build'
        run: 'hd'
        gpu: 'false'
        io: 'txt'
        grid: 'ug'
    """
    # Define default values for all possible arguments
    defaults = {
        'simulation': 'sedov',
        'hydro': 'rk2_gpu',
        'make': 'linux',
        'objdir': 'build',
        'run': 'hd',
        'gpu': 'false',
        'io': 'txt', 
        'grid': 'ug'
    }

    # Parse command-line arguments, overriding defaults
    for arg in sys.argv[1:]:
        # If user wants to know about configuration setting
        if arg.startswith('-') and (("help" in arg) or ("h" in arg)):
            print("Following are the setup default configurations:")
            print(defaults)
            print("Run as: ./setup +key=value")
            sys.exit(1)
        if arg.startswith('+') and '=' in arg:
            key, value = arg[1:].split('=', 1)
            if key in defaults:  # Only accept known keys
                defaults[key] = value
            else:
                print(f"Ignored +{key} = {value} as {key} is not recognized!")
    # Use correct subdirectory while using GPU
    #if defaults['gpu'] == "true":
    #    defaults['hydro'] = defaults['hydro'] + '_gpu'
    return defaults

def validate_paths(args):
    """
    Validate that required files and directories exist.
    Returns a tuple of validated paths and module implementations.
    """
    FilePath = "./"
    # Validate Makefile location
    site_dir = os.path.join(FilePath+'sites', args['make'])
    site_makefile = os.path.join(site_dir, 'Makefile')
    # if we are using gpu for our run
    if (args['gpu'] == "true"):
        site_makefile = os.path.join(site_dir, 'Makefile.gpu')
    if not os.path.exists(site_makefile):
        raise FileNotFoundError(f"Makefile not found: {site_makefile}")
    
    # Create object directory if needed
    objdir = FilePath+args['objdir']
    if os.path.exists(objdir):
        shutil.rmtree(objdir)
    os.makedirs(objdir)

    #Build module implementation map
    module_impls = {}
    
    for module, impl in args.items():
        if module not in ['make', 'objdir', 'run', 'gpu']:
            module_dir = os.path.join('src', module.upper())
            impl_dir = os.path.join(module_dir, f"{impl}")

            if not os.path.isdir(impl_dir):
                raise FileNotFoundError(
                    f"Implementation folder not found for {module}: {impl_dir}\n"
                )
            module_impls[module] = impl_dir
    
    return site_makefile, objdir, module_impls


def copy_src_files(objdir, module_impls, args):
    src_dir = "src"
    
    dest_dir = objdir
    for root, dirs, files in os.walk(src_dir):
        # Skip hidden directories
        dirs[:] = [d for d in dirs if not d.startswith('.')]
        # Create corresponding destination directory
        rel_path = os.path.relpath(root, src_dir)
        dest_path = os.path.join(objdir)

        # Copy all files from the COMMON directory
        if (rel_path == "COMMON"):
            for file in files:
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_path,file)
                shutil.copy(src_file, dest_file)
                print(f"Copied: {src_file}")


        #List of directories don't have any sub directory inside "src"
        NoSubDir = ["HD_h", "MHD_h"]
        Keys = ["run","run"]
        value = ["hd","mhd"]

        #List of directoires with sub-directory and implementation
        SubDirList = ["HYDRO", "IO", "GRID"]
        
        
        # Copy all the files from root directories
        if (root == src_dir):
            for file in files:
                
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_path, file)
                #print(src_file, dest_file)
                shutil.copy(src_file, dest_path)
                print(f"Copied: {src_file}")
        # Loop over these directories with no sub-directories
        # And copy over the files
        for kk in range(len(NoSubDir)):
            if ((args[Keys[kk]]) == value[kk]) and (rel_path == NoSubDir[kk]):
                for file in files:
                    
                    src_file = os.path.join(root, file)
                    dest_file = os.path.join(dest_path,file)
                    shutil.copy(src_file, dest_file)
                    print(f"Copied: {src_file}")
                    
        for ii in SubDirList:
            if (rel_path == ii.upper()):
                for file in files:
                    src_file = os.path.join(root, file)
                    dest_file = os.path.join(dest_path,file)
                    shutil.copy(src_file, dest_file)
                    print(f"Copied: {src_file}")
                for thisDir, subDir, files in os.walk(root):
                    if (thisDir == module_impls[ii.lower()]):
                        for file in files:
                            src_file = os.path.join(thisDir, file)
                            dest_file = os.path.join(dest_path,file)
                            shutil.copy(src_file, dest_file)
                            print(f"Copied: {src_file}")
                            
    # copy correct files from the problem                    
    for root, dirs, files in os.walk(src_dir):
        rel_path = os.path.relpath(root, src_dir)
        if (rel_path == "SIMULATION"): 
            for file in files: # copy files inside the src/SIMULATION
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_path,file)
                shutil.copy(src_file, dest_file)
                print(f"Copied: {src_file}")
        if (root == module_impls["simulation"]):
            for file in files: # copy files inside the src/SIMULATION/problem
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_path,file)
                shutil.copy(src_file, dest_file)
                print(f"Copied: {src_file}")      
    return

       

def main():
    try:
        # Parse arguments with defaults
        args = parse_arguments()
        print("Configuration settings:")
        for k, v in args.items():
            print(f"  {k}: {v}")

        # Validate paths and get module implementations
        site_makefile, objdir, module_impls = validate_paths(args)

        # Copy Makefile
        shutil.copy(site_makefile, os.path.join(objdir, 'Makefile'))
        print(f"Copied: {site_makefile} as Makefile")

        copy_src_files(objdir, module_impls, args)

        print("\nBuild configuration completed successfully!")
        print(f"Object directory: {objdir}")
        
    except Exception as e:
        print(f"\nError: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()


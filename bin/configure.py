import os
import sys
import shutil


def parse_arguments():
    # Default settings
    args = {
        'simulation': 'sedov',
        'hydro': 'rk2_gpu',
        'make': 'linux',
        'objdir': 'build',
        'run': 'hd',
        'ndim': '2',
        'gpu': 'false',
        'io': 'txt',
        'grid': 'ug'
    }

    for arg in sys.argv[1:]:
        if arg in ('-h', '--help'):
            print('Run as: ./setup +key=value')
            print('Options:', args)
            print('run: hd or mhd; ndim: 1, 2 or 3')
            sys.exit(0)
        if arg.startswith('+') and '=' in arg:
            key, value = arg[1:].split('=', 1)
            if key in args:
                args[key] = value
            else:
                raise ValueError(f'Unknown option: {key}')
        else:
            raise ValueError(f'Invalid argument: {arg}')

    if args['run'] not in ('hd', 'mhd'):
        raise ValueError('run must be hd or mhd')
    if args['ndim'] not in ('1', '2', '3'):
        raise ValueError('ndim must be 1, 2 or 3')
    if args['gpu'] not in ('true', 'false'):
        raise ValueError('gpu must be true or false')

    return args


def validate_paths(args):
    # Check inputs before clearing the build
    makefile = 'Makefile.gpu' if args['gpu'] == 'true' else 'Makefile'
    site_makefile = os.path.join('sites', args['make'], makefile)
    if not os.path.isfile(site_makefile):
        raise FileNotFoundError(f'Makefile not found: {site_makefile}')

    if not os.path.isfile('src/COMMON/param.h'):
        raise FileNotFoundError('Shared header not found: src/COMMON/param.h')

    module_impls = {}
    for module in ('hydro', 'io', 'grid', 'simulation'):
        impl_dir = os.path.join('src', module.upper(), args[module])
        if not os.path.isdir(impl_dir):
            raise FileNotFoundError(f'Implementation folder not found: {impl_dir}')
        module_impls[module] = impl_dir

    objdir = os.path.realpath(args['objdir'])
    repo = os.path.realpath('.')
    # Protect project files
    if os.path.islink(args['objdir']):
        raise ValueError('Build directory must not be a symlink')
    for protected in (repo, os.path.join(repo, 'src'),
                      os.path.join(repo, 'sites'), os.path.join(repo, 'bin'),
                      os.path.join(repo, '.git'), os.path.join(repo, 'notes'),
                      os.path.join(repo, 'par_files'), os.path.join(repo, 'test_files'),
                      os.path.join(repo, 'python_scripts'), os.path.join(repo, 'bash_scripts')):
        common = os.path.commonpath((objdir, protected))
        if common == objdir or (protected != repo and common == protected):
            raise ValueError(f'Unsafe build directory: {objdir}')

    return site_makefile, objdir, module_impls


def copy_files(source, objdir):
    for name in sorted(os.listdir(source)):
        src_file = os.path.join(source, name)
        if not name.startswith('.') and os.path.isfile(src_file):
            shutil.copy(src_file, os.path.join(objdir, name))
            print(f'Copied: {src_file}')


def copy_src_files(objdir, module_impls):
    copy_files('src', objdir)
    copy_files('src/COMMON', objdir)

    for module in ('hydro', 'io', 'grid', 'simulation'):
        copy_files(os.path.join('src', module.upper()), objdir)
        copy_files(module_impls[module], objdir)


def write_config(objdir, args):
    # Build settings
    lines = [
        f"#define NDIM {args['ndim']}"
    ]
    if args['run'] == 'mhd':
        lines.append('#define MHD')
    lines.append('')

    with open(os.path.join(objdir, 'config.h'), 'w') as config:
        config.write('\n'.join(lines))


def main():
    try:
        args = parse_arguments()
        print('Configuration settings:')
        for key, value in args.items():
            print(f'  {key}: {value}')

        site_makefile, objdir, module_impls = validate_paths(args)
        if os.path.exists(objdir):
            shutil.rmtree(objdir)
        os.makedirs(objdir)

        shutil.copy(site_makefile, os.path.join(objdir, 'Makefile'))
        copy_src_files(objdir, module_impls)
        write_config(objdir, args)

        print('\nBuild configuration completed successfully!')
        print(f'Object directory: {objdir}')
    except Exception as error:
        print(f'\nError: {error}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

import argparse
import sys
import os

from prism.DiagnoSR import DiagnoSR
from prism.a03_LookingGlass import LookingGlass
from prism.a04_ReGen import ReGen
from prism.paths import *
from prism.b01_utility import file_access
from prism.Interpreter import InterpreterInterface

# todo: Error handling for multiple genes in the same file in gene_databank


def get_function_branch(args):
    if args.data == 'gns':
        print("====[Browsing all available genes]====")
        for file in os.listdir(GENE_DATABANK):
            print(file)
    elif args.data == 'scr':
        print("====[Browsing all screened genes]====")
        for file in os.listdir(SCREEN_RESULTS):
            print(file)
    elif args.data == 'rpr':
        print(f"====[Browsing all repaired genes]====")
        for file in os.listdir(REGEN_CANDIDATES):
            print(file)
    else:
        print(f"Invalid PRISM data detected {args.data}")
        print(f"Choose from the following: [gns, scr, rpr]")

def access_filepath(args):
    try:
        file_access(args.branch, args.filename)
    except Exception as e:
        print(f"\nError during screening: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

def screen(args):
    print(f"[[Screening Initiated]]")
    output_filename = args.output if args.output is not None else None

    # file navigation diagnostic
    file_check, model, input_path, output_path = DiagnoSR(function_choice = "screen",
                          model_name=args.model,
                          input_filepath=args.input,
                          output_filename=output_filename).file_diagnostic()


    if file_check:
        try:
            LG = LookingGlass(model, genefile_route=input_path)
            print(f"\nScreening variant file: {input_path}")
            LG.predict_file(outfile_name=output_path)

        except Exception as e:
            print(f"\nError during screening: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    else:
        print(f"File check failed - ensure the following is true: {file_rules}")

def repair(args):
    print(f"[[Repair Initiated]]")
    # navigation
    output_filename = args.output if args.output is not None else None

    its = args.iterations if args.iterations is not None else 10
    cops = args.copies if args.copies is not None else 1

    # file navigation diagnostic
    file_check, model, input_path, output_path = DiagnoSR(function_choice = "repair",
                          model_name=args.model,
                          input_filepath=args.input,
                          output_filename=output_filename).file_diagnostic()

    if file_check:
        try:
            RG = ReGen(model,
                       pathogenic_route=input_path,
                       outfile_name=output_path,
                       iterations=its,
                       copies=cops)
            RG.repair()

        except Exception as e:
            print(f"\nError during repair: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    else:
        print(f"File check failed, ensure the following: {file_rules}")

def interact(args):
    print("Loading Interpreter...")

    verbose = False if not args.verbose else True
    model = args.model if args.model else "general"

    InterpreterInterface(args.data_type, args.dir_name, model=model, verbose=verbose).send_request()


def main():
    parser = argparse.ArgumentParser(
        prog='PRISM',  # name of program
        description="[Welcome to PRISM]"
                    "\nA derivative of GEM, specialized for variant pathogenicity screening and gene therapy target identification, complemented by a local AI-powered Interpretation module for result analysis and future scope propositions",
        formatter_class=argparse.RawDescriptionHelpFormatter,  # customizes help output
        epilog=  # text to display before argument
        """
        Parameters:
            note, parameters with "--" mean it's optional
        
            COMMAND LIST:
            [1] list
                - allows user to view available files in the following folders:
                
                'gns' -> what gene files we can screen or repair: from here we can select their names (minus the .fasta part) for input
                'scr' -> what genes have been screened
                'rpr' -> what genes have been put through ReGen
        
        
            [2 / 3] screen / repair 
                    - choice of screening for pathogenicity or repair
                
                input == input gene filepath  (no .fasta at the end)
                Note: if the input file is not in gene_databank you need to use its absolute path instead of input_filename
    
                --output == desired output filename if you don't want an auto-generated one
                
                --model clinicmod == Clinical Model
                --model discmod == Discriminator Model
    
                --iterations  == Regen Parameters / Run Options
                --copies
    
                output == output filename
            
            [4] access
                - allows user to look at any gene or result dataset
                datatype - ['gns', 'scr', 'rpr'] to select the folder
                dir_name - specific directory name to access -> use 'list' command to get list of all dir_names 
                note: for fasta files, do not type '.fasta'.
            
            [5] interpret
                - employs Interface class to use an AI model to interpret data on screened / repair results
                datatype - ['scr', 'rpr']
                dir_name
                --model_name -> will cause the AI to interpret results through these model's importance weights, blank will leave it to the top 15 features.
                --verbose -> see the prompt to be sent

        Command Examples:
            prism screen input_filename --output --model [clinicmod / discmod]
            prism repair input_filename --output --model [clinicmod / discmod] --copies x --iterations y
        """,
    )

    subparsers = parser.add_subparsers(dest='command', help='Command to execute')

    # ==============
    # list command - for seeing available files
    # ==============
    file_scouter = subparsers.add_parser(
        'list',
        help='Outputs all folder contents for future navigation ease'
    )

    file_scouter.add_argument(
        'data',
        type=str,
        help='input datatype -> gns: gene_databank, scr: screen_results, or rpr: regen_candidates'
    )

    # ==============
    # access command - for seeing file contents
    # ==============
    file_acc = subparsers.add_parser(
        'access',
        help='Outputs data in terminal'
    )

    file_acc.add_argument(
        'branch',
        type=str,
        help="choose from ['gns', 'scr', 'rpr']"
    )

    file_acc.add_argument(
        'filename',
        type=str,
        help='input filename'
    )


    # ==============
    # screen command
    # ==============
    screen_parser = subparsers.add_parser(
        'screen',
        help='Screen input variant for pathogenicity prediction',
        description='Screens input variant for pathogenicity prediction'
    )

    screen_parser.add_argument(
        'input',
        type=str,
        help='Input custom FASTA file w/ variants (provide filename inside gene_databank or full path)'
    )

    screen_parser.add_argument(
        '--model',
        type=str,
        default='clinicmod',
        choices=['clinicmod', 'discmod'],
        help='Model to use: "clinicmod" (lowest FNs), "discmod" (balanced pure discrimination)',
    )

    screen_parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Output filename prefix'
    )

    # ==============
    # repair command
    # ==============
    repair_parser = subparsers.add_parser(
        'repair',
        help='Repair input variant for pathogenicity reduction',
        description='Screen input variant for pathogenicity reduction & potential gene therapy target identification'
    )

    repair_parser.add_argument(
        'input',
        type=str,
        help='Input custom FASTA file w/ variants (provide filename inside gene_databank or full path)'
    )

    repair_parser.add_argument(
        '--model',
        type=str,
        default='clinicmod',
        choices=['clinicmod', 'discmod'],
        help='Model to use: "clinicmod" (lowest FNs), "discmod" (balanced pure discrimination)',
    )

    repair_parser.add_argument(
        '--iterations',
        type=int,
        default=10,
        help='Number of pathogenicity reduction iterations [default=10]'
    )

    repair_parser.add_argument(
        '--copies',
        type=int,
        default=1,
        help='Number of parallel sequence copies for pathogenicity reduction [default=1]'
    )

    repair_parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Output filename prefix'
    )

    # =================
    # interpret command - use local AI to interpret results
    # =================
    interpret = subparsers.add_parser(
        'interpret',
        help='Employs a local AI to interpret screen or repair results as well as generate future scope propositions'
    )

    interpret.add_argument(
        'data_type',
        type=str,
        help="choose from ['scr', 'rpr']"
    )

    interpret.add_argument(
        'dir_name',
        type=str,
        help='specific folder name to interpret'
    )

    interpret.add_argument(
        '--model',
        type=str,
        help='[clinicmod / discmod]'
    )

    interpret.add_argument(
        '--verbose',
        type=bool,
        help=''
    )


    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == 'screen':
        screen(args)
    elif args.command == 'repair':
        repair(args)
    elif args.command == 'list':
        get_function_branch(args)
    elif args.command == 'access':
        access_filepath(args)
    elif args.command == 'interpret':
        interact(args)


if __name__ == "__main__":
    main()

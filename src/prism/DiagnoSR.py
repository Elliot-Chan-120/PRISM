from pathlib import Path
from prism.paths import *
from prism.b01_utility import *

from prism.paths import REGEN_CANDIDATES, SCREEN_RESULTS


# file diagnostic tool for screen and repair functions
class DiagnoSR:
    def __init__(self, function_choice, model_name, input_filepath, output_filename=None):
        self.function_choice = function_choice
        self.input_filepath = input_filepath
        self.output_filename = output_filename
        self.model_name = model_name


    def file_diagnostic(self):
        print(f"===[[Initializing File Diagnostic]]===")
        self.validate_func()
        model = self.validate_model()
        input_path = self.validate_loc()
        output_path = self.validate_output()

        print("===[[Diagnostic Completed]]===")
        return True, model, input_path, output_path

    def validate_func(self):
        if self.function_choice == "screen" or "repair":
            return True
        else:
            raise ValueError("Invalid Function Provided")

    def validate_loc(self):
        """Validates filename in gene databank, or absolute path"""
        filename = f"{self.input_filepath}.fasta"
        print(f"validating input fasta file: {filename}")
        absolute_path = Path(filename)
        gene_databank_path = GENE_DATABANK / f"{filename}"

        if absolute_path.exists():
            print(f"file outside of gene_databank accepted | now loading {filename}")
            return str(absolute_path)
        elif gene_databank_path.exists():
            print(f"found in gene_databank | now loading {gene_databank_path}...")
            return str(gene_databank_path)
        else:
            # file not found
            print(f"Error: input filepath '{filename}' not found")
            print(f"Searched:")
            print(f" - {gene_databank_path}")
            print(f"{absolute_path}")
            print(f"\n ensure custom FASTA file is in correct folder location")
            return False

    def validate_output(self):
        """Validates output name, else automatically makes a new output filename based on function choice"""
        if self.output_filename is None:
            return f"{self.input_filepath}_{self.function_choice}ed"
        else:
            return self.output_filename

    def validate_model(self):
        """Validates model choice: discmod or clinicmod"""
        print(f"selecting model: {self.model_name}")
        if self.model_name == "clinicmod":
            return "ClinicalModel"
        elif self.model_name == "discmod":
            return "DiscriminatorModel"
        else:
            raise ValueError("Name must be one of: clinicmod / discmod")





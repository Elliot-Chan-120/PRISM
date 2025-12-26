from ollama import chat
from ollama import ChatResponse

from prism.paths import *
from prism.b01_utility import *

from prism.paths import REGEN_CANDIDATES, SCREEN_RESULTS

model_biollm = "gpt-oss"
# this is the model I downloaded off of ollama, but you are free to download a different one
# make sure to copy its EXACT NAME HERE

class InterpreterInterface:
    def __init__(self, data_type, dir_name, model, verbose=False):
        self.data_type = data_type

        if self.data_type == "rpr":
            self.request_type = 'repair_explanation'
            self.data_folder = REGEN_CANDIDATES
        elif self.data_type == "scr":
            self.request_type = 'screen_explanation'
            self.data_folder = SCREEN_RESULTS
        else:
            raise ValueError("Invalid mode detected")

        self.dir_name = dir_name
        self.dir_path = self.data_folder / dir_name

        self.model = model

        self.verbose = verbose


    def get_importances(self):
        importances = "\nHere are the model's top 20 feature importances among roughly 280. Their importances are recorded alongside their label\n"
        if "clinicmod" in self.model:
            importances += """
            FI_Mutation_composite: 0.0061
            Insertion_length: 0.0067
            2xScale_varState_delta__Intergenic_log2: 0.0070
            2xScale_varState_delta__5UTR_log2: 0.0073
            ap1_score: 0.0085
            tr_weighted_bploss: 0.0086
            FI_Mutation_propensity: 0.0095
            kozak_count: 0.0099
            tr_weighted_bpgain: 0.0138
            inr_count: 0.0163
            SNVs: 0.0178
            Global_score: 0.0223
            Stop_codons: 0.0243
            AA_Global_pct: 0.0251
            FI_General_Stability: 0.0294
            tr_contraction_score: 0.0338
            tr_expansion_score: 0.0630
            tr_composite_score: 0.0713
            tr_weighted_expansion_score: 0.0734
            tr_weighted_contraction_score: 0.0808
            tr_base_composite: 0.1113
            """

        elif "discmod" in self.model:
            importances += """
            2xScale_varState_delta__Intergenic_log2: 0.0067
            FI_relative_global: 0.0067
            relative_local: 0.0068
            ap1_score: 0.0070
            PTR_motif_count: 0.0071
            FI_Mutation_propensity: 0.0076
            creb_score: 0.0077
            1xScale_varState_delta__3UTR_log2: 0.0113
            inr_count: 0.0123
            Global_score: 0.0146
            FI_General_Stability: 0.0168
            tr_expansion_score: 0.0180
            AA_Global_pct: 0.0211
            Length_change: 0.0249
            Stop_codons: 0.0261
            tr_composite_score: 0.0518
            tr_weighted_contraction_score: 0.0608
            tr_weighted_expansion_score: 0.0616
            tr_contraction_score: 0.0623
            SNVs: 0.1066
            tr_base_composite: 0.1165
            """
        else:
            importances = "\nThese are the 15 features generally found to be of high importance across all model iterations. Generally, the most important features are at the bottom, but it may be important to note that they are amongst 280 total features.\n"
            importances += """
            Thermodynamic_stability
            FI_Mutation_propensity
            creb_score
            inr_count
            SNVs
            Global_score
            Stop_codons
            AA_Global_pct
            FI_General_Stability
            tr_expansion_score
            tr_composite_score
            tr_weighted_expansion_score
            tr_weighted_contraction_score
            tr_base_composite
            """

        return importances

    @staticmethod
    def search_folder(result_path, flag):
        files = os.listdir(result_path)
        for file in files:
            if flag in str(file):
                return result_path / file
        else:
            raise FileNotFoundError(f"{flag} file not in directory: re-run either scr or rpr")

    def assemble_prompt(self):
        prompt = """You are an expert computational biologist who has data on computationally generated gene therapy candidates. The program that acquired these candidates leverages ML models trained to recognize pathogenic and benign mutations by analyzing their biochemical / regulatory characteristic changes, what I call their Mutation Fingerprints. Explain, using the mutation fingerprints given, why the specific mutations produced their results (pathogenic / benign).\n"""
        # short feature category info
        prompt += """Feature categories include things like, regulatory motifs (CREB, inr, ap1, kozak), structural disruptions (stop codons, AA chain % similarity), repeat instability (tr_expansions, contractions .etc), biochemical features (thermodynamic stability) and feature interactions (FI). Note that percentages (% or pct) refer to composition, not probability."""

        prompt += self.get_importances()

        if self.request_type == 'repair_explanation':
            fingerprints = csv_parser(self.search_folder(self.dir_path, 'dataframe.csv'), verbose=False)
            mutation_log = txt_parser(self.search_folder(self.dir_path, 'results.txt'), verbose=False)
            unique = """\nYou're looking at results of a repair algorithm "regen" that iteratively generates benign computational candidates for pathogenic gene variants."""

            prompt += "\n===== MUTATION FINGERPRINTS =====\n"
            prompt += fingerprints
            prompt += "\n===== MUTATION LOG =====\n"
            prompt += mutation_log
            prompt += "\n"
            prompt += unique

        elif self.request_type == 'screen_explanation':
            fingerprints = csv_parser(self.search_folder(self.dir_path, 'screened_df.csv'))
            probabilities = csv_parser(self.search_folder(self.dir_path, 'screened.csv'))
            unique = """\nThese are the probabilities of the gene variant after being screened for pathogenicity."""

            prompt += "\n===== MUTATION FINGERPRINTS =====\n"
            prompt += fingerprints
            prompt += "\n===== PROBABILITIES =====\n"
            prompt += probabilities
            prompt += "\n"
            prompt += unique

        # reason with evidence
        prompt += """\nFor major changes in predicted benignity or pathogenicity, identify the most influential fingerprint features, explain how those changes plausibly contribute and note any uncertainty or conflicting signals."""
        # constraints
        prompt += """\nImportant constraints: Do not claim biological certainty. Treat all interpretations as hypotheses consistent with the provided features. If multiple explanations are plausible, present them as alternatives. Do not invent biologicla mechanisms not implied by the features"""

        # future scope
        prompt += """\nPropose: a computational follow-up, a wet-lab validation experiment, and one alternative repair strategy if the constraints prevent optimal changes"""


        if self.verbose:
            print(prompt)

        return prompt


    def send_request(self):
        print("Interacting with Interpreter...")
        try:
            response: ChatResponse = chat(
                model=model_biollm,
                messages=[{'role': 'user', 'content': self.assemble_prompt()}]
            )
            print(response['message']['content'])

        except Exception as e:
            print(f"Error while generating summary: {e}")


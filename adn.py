import tkinter as tk
from tkinter import filedialog, messagebox
import customtkinter
from PIL import ImageTk, Image
import random

# Customizing appearance using customtkinter

customtkinter.set_appearance_mode("dark")  # Modes: system (default), light, dark
customtkinter.set_default_color_theme("green")  # Themes: blue (default), dark-blue, green

class DNA:
    def __init__(self):
        self.sequence = ''

    def generate_random_sequence(self, length=100):
        self.sequence = ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(length))

    def load_sequence_from_file(self, filepath):
        with open(filepath, 'r') as file:
            self.sequence = ''.join(line.strip() for line in file if not line.startswith('>'))

    def validate_sequence(self):
        return set(self.sequence.upper()).issubset({'A', 'C', 'G', 'T'})

    def get_nucleotide_frequencies(self):
        return {nucleotide: self.sequence.count(nucleotide) for nucleotide in 'ACGT'}

    def transcribe_to_rna(self):
        return self.sequence.replace('T', 'U')

    def get_reverse_complement(self):
        COMPLEMENTS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(COMPLEMENTS[nucleotide] for nucleotide in reversed(self.sequence))

    def calculate_gc_content(self):
        return (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100

    def calculate_codon_frequencies(self):
        codons = [self.sequence[i:i+3] for i in range(0, len(self.sequence) - 2, 3)]
        codon_freq = {}
        for codon in codons:
            if codon in codon_freq:
                codon_freq[codon] += 1
            else:
                codon_freq[codon] = 1
        return codon_freq

    def translate_to_protein(self):
        CODON_TO_AMINOACID = {
            'AUG': 'Methionine', 'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine',
            'UUA': 'Leucine', 'UUG': 'Leucine', 'CUU': 'Leucine', 'CUC': 'Leucine',
            'CUA': 'Leucine', 'CUG': 'Leucine', 'AUU': 'Isoleucine', 'AUC': 'Isoleucine',
            'AUA': 'Isoleucine', 'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine',
            'GUG': 'Valine', 'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine',
            'UCG': 'Serine', 'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline',
            'CCG': 'Proline', 'ACU': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine',
            'ACG': 'Threonine', 'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine',
            'GCG': 'Alanine', 'UAU': 'Tyrosine', 'UAC': 'Tyrosine', 'CAU': 'Histidine',
            'CAC': 'Histidine', 'CAA': 'Glutamine', 'CAG': 'Glutamine', 'AAU': 'Asparagine',
            'AAC': 'Asparagine', 'AAA': 'Lysine', 'AAG': 'Lysine', 'GAU': 'Aspartic acid',
            'GAC': 'Aspartic acid', 'GAA': 'Glutamic acid', 'GAG': 'Glutamic acid',
            'UGU': 'Cysteine', 'UGC': 'Cysteine', 'UGG': 'Tryptophan', 'CGU': 'Arginine',
            'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine', 'AGU': 'Serine',
            'AGC': 'Serine', 'AGA': 'Arginine', 'AGG': 'Arginine', 'GGU': 'Glycine',
            'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine', 'UAA': 'Stop',
            'UAG': 'Stop', 'UGA': 'Stop'
        }
        protein = []
        for i in range(0, len(self.sequence) - 2, 3):
            codon = self.sequence[i:i+3]
            amino_acid = CODON_TO_AMINOACID.get(codon, '-')
            if amino_acid == 'Stop':
                break
            protein.append(amino_acid)
        return protein

    def calculate_protein_mass(self):
        AMINOACID_WEIGHTS = {
            'Methionine': 149, 'Phenylalanine': 165, 'Leucine': 131, 'Isoleucine': 131,
            'Valine': 117, 'Serine': 105, 'Proline': 115, 'Threonine': 119, 'Alanine': 89,
            'Tyrosine': 181, 'Histidine': 155, 'Glutamine': 146, 'Asparagine': 132,
            'Lysine': 146, 'Aspartic acid': 133, 'Glutamic acid': 147, 'Cysteine': 121,
            'Tryptophan': 204, 'Arginine': 174, 'Glycine': 75
        }
        protein_mass = sum(AMINOACID_WEIGHTS.get(amino_acid, 0) for amino_acid in self.translate_to_protein() if amino_acid != '-')
        nucleotide_freq = self.get_nucleotide_frequencies()
        return protein_mass, nucleotide_freq

class DNAApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Analyseur de DNA")
        self.geometry("600x400")
        self.dna = DNA()
        
        img1 = ImageTk.PhotoImage(Image.open("./pattern.png"))
        l1 = customtkinter.CTkLabel(master=self, image=img1)
        l1.pack()
        self.create_widgets()

    

    def create_widgets(self):
        # Custom frame
        frame = customtkinter.CTkFrame(master=self, width=320, height=500, corner_radius=0)
        frame.place(relx=0.5, rely=0.4, anchor=tk.CENTER)

        # Custom buttons
        button1 = customtkinter.CTkButton(master=frame, width=220, text="Générer ADN aléatoire", command=self.generate_dna, corner_radius=6)
        button1.place(x=50, y=30)

        button2 = customtkinter.CTkButton(master=frame, width=220, text="Charger ADN depuis fichier", command=self.load_dna, corner_radius=6)
        button2.place(x=50, y=80)

        button3 = customtkinter.CTkButton(master=frame, width=220, text="Valider ADN", command=self.validate_dna, corner_radius=6)
        button3.place(x=50, y=130)

        button4 = customtkinter.CTkButton(master=frame, width=220, text="Transcrire en ARN", command=self.transcribe_rna, corner_radius=6)
        button4.place(x=50, y=180)

        button5 = customtkinter.CTkButton(master=frame, width=220, text="Calculer complément inverse", command=self.calculate_reverse_complement, corner_radius=6)
        button5.place(x=50, y=230)

        button6 = customtkinter.CTkButton(master=frame, width=220, text="Calculer le taux de GC", command=self.calculate_gc_content, corner_radius=6)
        button6.place(x=50, y=280)

        button7 = customtkinter.CTkButton(master=frame, width=220, text="Calculer les fréquences de codons", command=self.calculate_codon_frequencies, corner_radius=6)
        button7.place(x=50, y=330)

        button8 = customtkinter.CTkButton(master=frame, width=220, text="Traduire en protéine", command=self.translate_to_protein, corner_radius=6)
        button8.place(x=50, y=380)

        button9 = customtkinter.CTkButton(master=frame, width=220, text="Calculer la masse protéique", command=self.calculate_protein_mass_display, corner_radius=6)
        button9.place(x=50, y=430)

        # Result text display
        self.result_text = tk.Text(self, height=8, width=40, bg="black", fg="white", font=("Arial", 12))
        self.result_text.place(relx=0.5, rely=0.85, anchor=tk.CENTER)


    def generate_dna(self):
        self.dna.generate_random_sequence(100)
        self.update_result(f"ADN généré: {self.dna.sequence}")

    def load_dna(self):
        filepath = filedialog.askopenfilename()
        if filepath:
            self.dna.load_sequence_from_file(filepath)
            self.update_result(f"ADN chargé: {self.dna.sequence}")

    def validate_dna(self):
        if self.dna.validate_sequence():
            messagebox.showinfo("Validation", "La séquence ADN est valide.")
        else:
            messagebox.showerror("Validation", "La séquence ADN est invalide.")

    def transcribe_rna(self):
        rna = self.dna.transcribe_to_rna()
        self.update_result(f"ARN transcrit: {rna}")

    def calculate_reverse_complement(self):
        reverse_complement = self.dna.get_reverse_complement()
        self.update_result(f"Complément inverse: {reverse_complement}")

    def calculate_gc_content(self):
        gc_content = self.dna.calculate_gc_content()
        self.update_result(f"Taux de GC: {gc_content}%")

    def calculate_codon_frequencies(self):
        codon_freq = self.dna.calculate_codon_frequencies()
        self.update_result(f"Fréquences de codons: {codon_freq}")

    def translate_to_protein(self):
        protein = self.dna.translate_to_protein()
        self.update_result(f"Protéine traduite: {' '.join(protein)}")

    def calculate_protein_mass_display(self):
        protein_mass, nucleotide_freq = self.dna.calculate_protein_mass()
        result = f"Masse protéique: {protein_mass:.2f} g/mol\n"
        result += f"Fréquence de chaque nucléotide:\n"
        for nucleotide, freq in nucleotide_freq.items():
            result += f"{nucleotide}: {freq}\n"
        self.update_result(result)

    def update_result(self, message):
        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, message)


if __name__ == "__main__":
    app = DNAApp()
    app.mainloop()

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from grocharmm.gro_to_crd import CRDeditor
from grocharmm.gro_to_psf import PSFeditor
from grocharmm.utils import extract_segs_crd, DEFAULT_SEGMENT_RULES

class EditorGUI:
    def __init__(self, root):
        self.root = root
        root.title("CHARMM File Editor")

        self.mode = tk.StringVar(value="crd")
        self.input_path = tk.StringVar()
        self.output_path = tk.StringVar()

        # File Selection
        file_frame = ttk.LabelFrame(root, text="File Selection")
        file_frame.grid(row=0, column=0, columnspan=3, padx=10, pady=10, sticky="ew")

        ttk.Label(file_frame, text="Input file").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        ttk.Entry(file_frame, textvariable=self.input_path, width=40).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(file_frame, text="Browse", command=self.select_input).grid(row=0, column=2, padx=5, pady=5)

        ttk.Label(file_frame, text="Output file").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        ttk.Entry(file_frame, textvariable=self.output_path, width=40).grid(row=1, column=1, padx=5, pady=5)
        ttk.Button(file_frame, text="Browse", command=self.select_output).grid(row=1, column=2, padx=5, pady=5)

        ttk.Button(file_frame, text="Preview Segments", command=self.preview_segments).grid(row=2, column=1, columnspan=2, pady=5)

        # Mode
        mode_frame = ttk.LabelFrame(root, text="Mode")
        mode_frame.grid(row=1, column=0, columnspan=3, padx=10, pady=(0,10), sticky="ew")
        ttk.Radiobutton(mode_frame, text="CRD", variable=self.mode, value="crd").grid(row=0, column=0, padx=5, pady=5)
        ttk.Radiobutton(mode_frame, text="PSF", variable=self.mode, value="psf").grid(row=0, column=1, padx=5, pady=5)

        # Segment Rules
        segment_frame = ttk.LabelFrame(root, text="Custom Segment Rules (resname:segment)")
        segment_frame.grid(row=2, column=0, columnspan=3, padx=10, pady=(0, 10), sticky="ew")
        self.segment_rules_entry = tk.Text(segment_frame, height=6, width=60)

        default_text = "\n".join(f"{k}:{v}" for k, v in sorted(DEFAULT_SEGMENT_RULES.items()))
        self.segment_rules_entry.insert("1.0", default_text)
        self.segment_rules_entry.grid(row=0, column=0, padx=5, pady=5)

        # Options
        self.replace_resnames = tk.BooleanVar()
        tk.Checkbutton(root, text="Replace resnames (CRD only)", variable=self.replace_resnames).grid(row=3, columnspan=2)

        # Run
        tk.Button(root, text="Run", command=self.run_editor).grid(row=4, column=1)

    def select_input(self):
        path = filedialog.askopenfilename()
        if path:
            self.input_path.set(path)

    def select_output(self):
        ext = ".crd" if self.mode.get() == "crd" else ".psf"
        path = filedialog.asksaveasfilename(defaultextension=ext)
        if path:
            self.output_path.set(path)

    def parse_segment_rules(self):
        raw_text = self.segment_rules_entry.get("1.0", "end").strip()
        rules = {}
        for line in raw_text.splitlines():
            if ":" in line:
                resname, segname = line.split(":", 1)
                rules[resname.strip()] = segname.strip()
        return rules

    def preview_segments(self):
        input_path = self.input_path.get()
        if not input_path:
            messagebox.showerror("No input file", "Please select an input CRD file first.")
            return

        try:
            segments = extract_segs_crd(input_path)
            if not segments:
                messagebox.showinfo("Segments Found", "No segments found in the file.")
                return

            # Get current entries
            current_text = self.segment_rules_entry.get("1.0", "end").strip()
            current_rules = {}
            for line in current_text.splitlines():
                if ":" in line:
                    resname, segname = line.split(":", 1)
                    current_rules[resname.strip()] = segname.strip()

            # Fill in missing with defaults or '?'
            for seg in segments:
                if seg not in current_rules:
                    current_rules[seg] = DEFAULT_SEGMENT_RULES.get(seg, "?")

            updated_text = "\n".join(f"{k}:{v}" for k, v in sorted(current_rules.items()))
            self.segment_rules_entry.delete("1.0", "end")
            self.segment_rules_entry.insert("1.0", updated_text)

            messagebox.showinfo("Segments Found", f"Found segments:\n\n{', '.join(sorted(segments))}")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to extract segments:\n{e}")

    def run_editor(self):
        input_path = self.input_path.get()
        output_path = self.output_path.get()

        if not input_path or not output_path:
            messagebox.showerror("Missing file paths", "Please select both input and output files.")
            return

        segment_rules = self.parse_segment_rules()

        try:
            if self.mode.get() == "crd":
                editor = CRDeditor(input_path)
                editor.read()
                if self.replace_resnames.get():
                    editor.replace_resnames()
                editor.update_segments(segment_rules=segment_rules)
                editor.update_segment_numbers()
                editor.write(output_path)
            else:
                editor = PSFeditor(input_path)
                editor.read_lines()
                editor.update_psf_segments(segment_rules=segment_rules)
                editor.update_psf_resids()
                editor.write_lines(output_path)

            messagebox.showinfo("Success", f"File saved to:\n{output_path}")

        except Exception as e:
            messagebox.showerror("Error", f"Something went wrong:\n{e}")

def main():
    root = tk.Tk()
    gui = EditorGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()

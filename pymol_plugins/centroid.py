from pymol import cmd
import numpy as np

def render_centroid(selection="sele"):
    """
    Computes the centroid of the selected residues and renders a pseudoatom at the centroid.
    
    Args:
        selection (str): PyMOL selection string specifying the residues of interest.
    """
    # Get atomic coordinates of selected residues
    model = cmd.get_model(selection)
    
    # Extract coordinates as a numpy array
    coords = np.array([atom.coord for atom in model.atom])

    # Check if any atoms were found
    if coords.shape[0] == 0:
        print("Error: No atoms found in the selection.")
        return
    
    # Compute centroid
    centroid = coords.mean(axis=0)

    # Ensure the centroid is a tuple of floats
    centroid_tuple = tuple(map(float, centroid))  # Convert numpy array to Python tuple

    # Create a pseudoatom at the centroid
    cmd.pseudoatom("centroid", pos=centroid_tuple, color="red", label="Centroid")

    # Show the centroid as a sphere
    cmd.show("spheres", "centroid")
    cmd.set("sphere_scale", 1.0, "centroid")

    print(f"Centroid coordinates: {centroid_tuple}")

# Example Usage: Run in PyMOL
render_centroid("chain A and resi 50+52+60")

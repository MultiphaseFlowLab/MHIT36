import vtk
import os

# PARAMETERS
data_dir = "/leonardo_scratch/large/userexternal/aroccon0/VANDYKE/post/paraview_hit/output"
file_pattern = "OUTPAR_*.vtk"
iso_value = 0.5
output_dir = "stl"
array_name = "phi"

os.makedirs(output_dir, exist_ok=True)
files = sorted(f for f in os.listdir(data_dir) if f.endswith(".vtk"))

for fname in files:
    step = fname.split("_")[-1].replace(".vtk", "")
    filepath = os.path.join(data_dir, fname)

    print(f"Processing step {step}")

    # -----------------------------
    # READ VTK FILE
    # -----------------------------
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(filepath)
    reader.Update()

    data = reader.GetOutput()

    # -----------------------------
    # SELECT SCALAR ARRAY
    # -----------------------------
    if data.GetPointData().HasArray(array_name):
        data.GetPointData().SetActiveScalars(array_name)
        print("  phi on POINTS")
    elif data.GetCellData().HasArray(array_name):
        data.GetCellData().SetActiveScalars(array_name)
        print("  phi on CELLS")
    else:
        raise RuntimeError(f"{array_name} not found in {fname}")

    # -----------------------------
    # ISO-SURFACE (ParaView Contour)
    # -----------------------------
    contour = vtk.vtkContourFilter()
    contour.SetInputData(data)
    contour.SetValue(0, iso_value)
    contour.Update()

    poly = contour.GetOutput()
    n_cells = poly.GetNumberOfCells()
    print(f"  contour cells: {n_cells}")

    if n_cells == 0:
        print("  ⚠️ empty iso-surface, skipping")
        continue

    # -----------------------------
    # WRITE STL
    # -----------------------------
    outname = os.path.join(output_dir, f"phi_iso_{step}.stl")

    writer = vtk.vtkSTLWriter()
    writer.SetFileName(outname)
    writer.SetInputData(poly)
    writer.SetFileTypeToBinary()
    writer.Write()

print("✅ All STL files written.")

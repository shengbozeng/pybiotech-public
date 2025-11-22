
from pathlib import Path
from pycorelibs.network.requests import HTTPMethod, fetch_url
from pybiotech.loaders.sdf_loader import SDFLoader




if __name__ == "__main__":
        
    sdf_file = Path("C:\\迅雷下载\\69553331_69580545.sdf")
    # sdf_file = Path("D:\\all_data_from_pubchem.2024-12-11_20-39-43.sdf")
    # sdf_file = Path("C:\\projects\\MolecularViewer\\69553331_69580545.sdf")
    mol_list,t,a, = SDFLoader.readDataFromFile(str(sdf_file))

    for i, entry in enumerate(mol_list):
        if i>10:
            break
        print(f"[{i}] {entry.GetPropsAsDict()}")
        print('*'*80)
        pass
from typing import List,Dict
from rich.progress import Progress
from pybiotech.loaders.nih.pubchem.online.compound import (
    get_compound,
    get_similarity_compound,
)
from pybiotech.type.nih.pubchem import (
    ALNPCompound,
    EInputType,
    EOperationType,
    EOutputType,
)


test_data = [
    100203,
    155245,
    11226684,
    9868928,
    10141,
    3806,
    8530,
    69457,
    8342,
    137321687,
    16861,
    74944,
]


def test_get_compound(cid_list: List[str]) -> Dict[str, ALNPCompound]:
    with Progress() as progress:
        task = progress.add_task(
            "[cyan]Get PubChem online data...", total=1
        )  # 先放1，回调里会动态更新 total

        def cb(done: int, total: int, msg: str):
            # 动态更新 total，显示 done/total 和百分比
            progress.update(
                task,
                total=total,
                completed=done,
                description=f"[cyan]{msg}  ({done}/{total})",
            )

        data = get_compound(
            cid_list=cid_list, include_conformer=True, progress_callback=cb
        )
    return data


def test_get_similarity_compound(cid_list: List[str]) -> Dict[str, ALNPCompound]:
    with Progress() as progress:
        task = progress.add_task(
            "[cyan]Get PubChem online data...", total=1
        )  # 先放1，回调里会动态更新 total

        def cb(done: int, total: int, msg: str):
            # 动态更新 total，显示 done/total 和百分比
            progress.update(
                task,
                total=total,
                completed=done,
                description=f"[cyan]{msg}  ({done}/{total})",
            )

        data = get_similarity_compound(
            input=EInputType.CID,
            value=cid_list[-3],
            operation=EOperationType.CIDS,
            output=EOutputType.TXT,
            include_conformer=True,
            progress_callback=cb,
        )
    return data


if __name__ == "__main__":

    cid_list = [str(cid) for cid in test_data]
    # data = test_get_compound(cid_list=cid_list)

    data = test_get_similarity_compound(cid_list=cid_list)
    for cid, compound in data.items():
        print(f"CID: {cid}")
        print("Compound Data:\n", compound)
        print('*' * 50)
    print("✅ 全部完成")

    # for cid in cid_list:
    #     print(f"CID: {cid}")
    #     print(get_similarity_by_smiles(input=EInputType.CID, value=cid, operation=EOperationType.CIDS, output=EOutputType.JSON))

    # smiles = "CCCCCC1C(C(OC(=O)C(C(OC1=O)C)NC(=O)C2=C(C(=CC=C2)NC=O)O)C)OC(=O)CC(C)C"

    # print(get_similarity(input=EInputType.CID, value=cid_list[-1], operation=EOperationType.CIDS, output=EOutputType.TXT))

    pass

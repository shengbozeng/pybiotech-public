# -*- coding: utf-8 -*-
''' ************************************************************ 
### Author: Zeng Shengbo shengbo.zeng@ailingues.com
### Date: 2025-06-26 14:55:40
### LastEditors: Zeng Shengbo shengbo.zeng@ailingues.com
### LastEditTime: 2025-07-01 09:34:31
### FilePath: //ailingues//test//test_uniprot_loader.py
### Description: 
### 
### Copyright (c) 2025 by AI Lingues, All Rights Reserved. 
********************************************************** '''
from pycorelibs.network.requests import HTTPMethod, fetch_url
from pybiotech.loaders.uniprot_loader import UniProtLoader

# xml_file = Path("E:\\迅雷下载\\uniprot_sprot.xml\\uniprot_sprot.xml")
# loader = UniProtLoader(xml_file)

# for i, entry in enumerate(loader.iterate_entries()):
#     if i>100:
#         break
#     print(f"[{i}] {entry.accession[0]} - {entry.name[0]} - Protein:\t{entry.protein.recommended_name.full_name}")
#     if i % 500 == int(i/500):
#         pass
#     pass

url = "https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)%20AND%20(organism_id:9606)&format=xml"
if __name__ == "__main__":
    # asyncio.run(main())
    result = fetch_url(
        url=url,
        method=HTTPMethod.GET,
        headers={"Content-Type": "application/xml"},
        # json_data={"message": "hello"},
        # proxies={"http": "http://127.0.0.1:8888"},
        max_retries=2,
    )

    if result["success"]:
        print("请求成功：", result["content"])
        loader = UniProtLoader(result["content"])

        for i, entry in enumerate(loader.iterate_entries()):
            print(f"[{i}] {entry.accession[0]} - {entry.name[0]} - Protein:\t{entry.protein.recommended_name.full_name}")
            pass
        pass
    else:
        print("请求失败：", result["error"])
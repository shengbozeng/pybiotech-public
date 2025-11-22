# -*- coding: utf-8 -*-
''' ************************************************************ 
### Author: Zeng Shengbo shengbo.zeng@ailingues.com
### Date: 2025-06-26 20:53:20
### LastEditors: Zeng Shengbo shengbo.zeng@ailingues.com
### LastEditTime: 2025-07-02 09:12:04
### FilePath: //ailingues//test//test_query_field_loader.py
### Description: 
### 
### Copyright (c) 2025 by AI Lingues, All Rights Reserved. 
********************************************************** '''

import argparse
from pybiotech.loaders.uniprot_query_field_loader import UniProtQueryFieldLoader
from pycorelibs.network.requests import HTTPMethod, fetch_url


def main():
    import os,sys
    from pathlib import Path
    fields_json_path =  Path(os.path.dirname(__file__)) / 'UniProtQueryFieldsDefinition.json'
    url = 'https://rest.uniprot.org/configure/uniprotkb/search-fields'
    result = fetch_url(
        url=url,
        method=HTTPMethod.GET,
        headers={"Content-Type": "application/json"},
        # json_data={"message": "hello"},
        # proxies={"http": "http://127.0.0.1:8888"},
        max_retries=2,
    )

    if result["success"]:
        print("请求成功：", result["content"])
    else:
        print("请求失败：", result["error"])
        
    parser = argparse.ArgumentParser(description="构造 UniProt 查询语句")
    parser.add_argument("--json", required=False,default=result["content"], help="字段定义 JSON 路径")
    # parser.add_argument("--term", required=False, default='sec_acc', help="字段 term，如 gene 或 accession")
    # parser.add_argument("--value", required=False, default='B2R5V1', help="查询值，如 P53 或 P12345")

    args = parser.parse_args()
    term_value_dict={
        "reviewed":'true',
        "gene":"PSMB5"
    }
    fields = UniProtQueryFieldLoader.load_query_fields(args.json)   
    for k,v in term_value_dict.items(): 
        term_field = UniProtQueryFieldLoader.find_field_by_term(fields, k)
        if not term_field:
            print(f"[错误] 未找到字段 term = '{k}'")
            sys.exit(1)
            
        query = UniProtQueryFieldLoader.build_query(term_field, v)
        print(query)  # 输出：gene:P53


if __name__ == "__main__":
    main()

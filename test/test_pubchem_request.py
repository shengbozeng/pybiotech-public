# -*- coding: utf-8 -*-
''' ************************************************************ 
### Author: Zeng Shengbo shengbo.zeng@ailingues.com
### Date: 2025-07-17 16:11:48
### LastEditors: Zeng Shengbo shengbo.zeng@ailingues.com
### LastEditTime: 2025-07-17 16:32:18
### FilePath: //ailinguesbiotechlib//test//test_pubchem_request.py
### Description: 
### 
### Copyright (c) 2025 by AI Lingues, All Rights Reserved. 
********************************************************** '''
from pybiotech.classes.nih.https.pubchem.input_spec import InputSpec
from pybiotech.classes.nih.https.pubchem.opera_spec import OperationSpec
from pybiotech.classes.nih.https.pubchem.output_spec import OutputSpec
from pybiotech.classes.nih.https.pubchem.query_options import QueryOptions
from pybiotech.classes.nih.https.pubchem.request import PubChemRequest

def demo_success():
    print("【Demo: 正确请求】")
    # 例1：查询 aspirin 分子的分子式、InChIKey，JSON输出
    try:
        input_spec = InputSpec(
            domain="compound",
            namespace="name",
            identifiers="aspirin"
        )
        operation_spec = OperationSpec(
            operation="property",
            tags=["MolecularFormula", "InChIKey"],
            domain="compound"  # 推荐补上，校验更严
        )
        output_spec = OutputSpec(output_format="JSON")
        query_options = QueryOptions(list_return="flat", MaxRecords=5)
        # 设置上下文进行参数合法性校验
        query_options.set_context(operation_spec.operation, output_spec.output_format)

        request = PubChemRequest(input_spec, operation_spec, output_spec, query_options)
        url = request.build_url()
        headers = request.build_headers()
        print("API URL:", url)
        print("Headers:", headers)
        print("所有参数校验通过，准备请求PubChem。")
    except Exception as e:
        print("参数错误：", e)

def demo_error_mutual_exclusion():
    print("\n【Demo: 错误互斥校验】")
    # 例2：错误用法——callback和非JSONP混用，应该报错
    try:
        input_spec = InputSpec(
            domain="compound",
            namespace="cid",
            identifiers="2244"
        )
        operation_spec = OperationSpec(
            operation="synonyms",
            tags=["MolecularFormula", "InChIKey"],  # 明确 tags 非空
            domain="compound"
        )
        output_spec = OutputSpec(output_format="JSON")  # 注意不是JSONP
        # 错误地添加 callback
        query_options = QueryOptions(callback="some_callback_name")
        query_options.set_context(operation_spec.operation, output_spec.output_format)

        request = PubChemRequest(input_spec, operation_spec, output_spec, query_options)
        url = request.build_url()  # 这里应该会抛出参数冲突异常
        print("API URL:", url)
        print("Headers:", request.build_headers())
    except Exception as e:
        print("捕获到参数冲突：", e)

def demo_error_missing_required():
    print("\n【Demo: 缺失必需参数校验】")
    # 例3：property操作但没指定tags
    try:
        input_spec = InputSpec(
            domain="compound",
            namespace="cid",
            identifiers="2244"
        )
        operation_spec = OperationSpec(
            operation="property",  # property 必须带tags
            domain="compound"
        )
        output_spec = OutputSpec(output_format="JSON")
        query_options = QueryOptions()
        query_options.set_context(operation_spec.operation, output_spec.output_format)

        request = PubChemRequest(input_spec, operation_spec, output_spec, query_options)
        url = request.build_url()  # 这里应该会报property必须有tags
        print("API URL:", url)
    except Exception as e:
        print("捕获到缺失必需参数：", e)

if __name__ == "__main__":
    demo_success()
    demo_error_mutual_exclusion()
    demo_error_missing_required()

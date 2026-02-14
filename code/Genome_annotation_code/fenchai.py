def split_genbank(input_file, output_prefix):
    with open(input_file, 'r') as f:
        content = f.read()
    
    # 按 "//" 分割记录（注意：可能有空行，需处理）
    records = [r.strip() for r in content.split('//') if r.strip()]
    
    for i, record in enumerate(records, 1):
        # 每个记录以 "//" 结尾，但我们需要去掉它（因为 split 已经分割了）
        full_record = record + '//\n'  # 恢复结尾的 "//"
        output_file = f"{output_prefix}_record_{i}.gbff"
        with open(output_file, 'w') as f:
            f.write(full_record)
        print(f"Saved: {output_file}")

# 使用示例
split_genbank('benzhouyidingzuowanzhushi_GLOBAL_multi-GenBank.gbff', 'split_record')

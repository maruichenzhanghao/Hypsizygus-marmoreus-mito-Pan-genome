from Bio import SeqIO
import argparse

def fasta_length(file_path):
    """
    读取FASTA文件并打印每条序列的ID及其长度。
    
    :param file_path: FASTA文件路径
    """
    try:
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                print(f"Sequence ID: {record.id}, Length: {len(record.seq)}")
    except FileNotFoundError:
        print(f"The file {file_path} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # 创建解析器对象
    parser = argparse.ArgumentParser(description="Calculate the length of each sequence in a FASTA file.")
    
    # 添加命令行参数
    parser.add_argument("fasta_file", type=str, help="Path to the FASTA file")
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 使用传入的参数调用函数
    fasta_length(args.fasta_file)

singularity exec /home/maxinxin/software/01.singlity/cactus.sif \
cactus-pangenome \
  ./jobstore \          # 工作目录
  ./mito_samples.txt \  # 样本文件
  --outDir mito_results \      # 输出目录
  --outName mito_pg \          # 输出前缀
  --reference f2 \# 参考样本名
  
  # 线粒体专用核心参数
  --noSplit \                  # 保持序列完整
  --maxLen 200000 \             # 覆盖最大长度
  --clip 0 \                   # 禁用末端修剪
  
  # 输出格式（按您需求）
  --vcf \                      # VCF变异文件
  --gfa full \                 # 完整GFA图（含路径）
  --gbz \                      # GBZ压缩格式
  --xg \                       # vg兼容格式
  --odgi \                     # 实验性ODGI输出
  --viz \                      # 自动可视化
  
  # 资源控制（示例值，需调整）
  --mgCores 16 \               # minigraph核心数
  --mgMemory 16G \             # minigraph内存（线粒体需较少）
  --indexCores 8 \             # 索引构建核心
  --indexMemory 8G \           # 索引内存
  
  # 高级优化
  --vcfwave \                  # 改进变异检测
  --giraffe \                  # 生成Giraffe索引
  --lastTrain \                # 优化比对参数

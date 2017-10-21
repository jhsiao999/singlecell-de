library(Biobase)
library(ashbun)
library(singleCellRNASeqMouseZeiselBrain)

# --- Save processed/cleaned raw data
|-- /project2/gilad/joycehsiao/singlecell-de
|    |-- data/
|        |-- mousezeiselbrain.rds
|        |-- humantungipsc.rds
|        |-- mousekleinesc.rds
|        |-- gtex001.rds: liver samples with thinning paramter .001
|        |-- gtex01.rds: liver samples with thinning paramter .01
|        |-- gtex.rds: liver samples
|   |-- output_eval/
|       |-- mousezeiselbrain/
|           |-- mousezeiselbrain.allgenes.pi05.n50.bignormal.output.rds
|           |-- mousezeiselbrain.allgenes.pi09.n50.bignormal.output.rds
|           |-- mousezeiselbrain.pergene.pi05.n50.bignormal.output.rds
|           |-- mousezeiselbrain.pergene.pi09.n50.bignormal.output.rds
|       |-- humantungipsc/
|           |-- humantungipsc.allgenes.pi05.n50.bignormal.output.rds
|           |-- humantungipsc.allgenes.pi09.n50.bignormal.output.rds
|           |-- humantungipsc.pergene.pi05.n50.bignormal.output.rds
|           |-- humantungipsc.pergene.pi09.n50.bignormal.output.rds
|       |-- mousekleinesc/
|           |-- mousekleinesc.allgenes.pi05.n50.bignormal.output.rds
|           |-- mousekleinesc.allgenes.pi09.n50.bignormal.output.rds
|           |-- mousekleinesc.pergene.pi05.n50.bignormal.output.rds
|           |-- mousekleinesc.pergene.pi09.n50.bignormal.output.rds
|       |-- gtex001/
|           |-- gtex001.allgenes.pi05.n50.bignormal.output.rds
|           |-- gtex001.allgenes.pi09.n50.bignormal.output.rds
|           |-- gtex001.pergene.pi05.n50.bignormal.output.rds
|           |-- gtex001.pergene.pi09.n50.bignormal.output.rds
|   |-- output_rocavg/
|       |-- mousezeiselbrain/
|           |-- mousezeiselbrain.allgenes.pi05.n50.bignormal.rocavg.rds
|           |-- mousezeiselbrain.allgenes.pi09.n50.bignormal.rocavg.rds
|           |-- mousezeiselbrain.pergene.pi05.n50.bignormal.rocavg.rds
|           |-- mousezeiselbrain.pergene.pi09.n50.bignormal.rocavg.rds
|       |-- humantungipsc/
|           |-- humantungipsc.allgenes.pi05.n50.bignormal.rocavg.rds
|           |-- humantungipsc.allgenes.pi09.n50.bignormal.rocavg.rds
|           |-- humantungipsc.pergene.pi05.n50.bignormal.rocavg.rds
|           |-- humantungipsc.pergene.pi09.n50.bignormal.rocavg.rds
|       |-- mousekleinesc/
|           |-- mousekleinesc.allgenes.pi05.n50.bignormal.rocavg.rds
|           |-- mousekleinesc.allgenes.pi09.n50.bignormal.rocavg.rds
|           |-- mousekleinesc.pergene.pi05.n50.bignormal.rocavg.rds
|           |-- mousekleinesc.pergene.pi09.n50.bignormal.rocavg.rds
|       |-- gtex.001/
|           |-- gtex.001.allgenes.pi05.n50.bignormal.rocavg.rds
|           |-- gtex.001.allgenes.pi09.n50.bignormal.rocavg.rds
|           |-- gtex.001.pergene.pi05.n50.bignormal.rocavg.rds
|           |-- gtex.001.pergene.pi09.n50.bignormal.rocavg.rds



# ---- Make simulated data and save in
|-- /scratch/midway2/joycehsiao/singlecell-de
|   |-- data_n50/
|       |-- filter01/
|           |-- mousezeiselbrain.n50.filter01.rds
|           |-- humantungipsc.n50.filter01.rds
|           |-- mousekleinesc.n50.filter01.rds
|           |-- gtex001.n50.filter01.rds
|   |-- data_n50_filter01_simulated/
|       |-- mousezeiselbrain/
|           |-- mousezeiselbrain.allgenes.filter01.pi05.n50.bignormal.data.rds
|           |-- mousezeiselbrain.allgenes.filter01.pi09.n50.bignormal.data.rds
|           |-- mousezeiselbrain.pergene.filter01.pi05.n50.bignormal.data.rds
|           |-- mousezeiselbrain.pergene.filter01.pi09.n50.bignormal.data.rds
|       |-- humantungipsc/
|           |-- humantungipsc.allgenes.filter01.pi05.n50.bignormal.data.rds
|           |-- humantungipsc.allgenes.filter01.pi09.n50.bignormal.data.rds
|           |-- humantungipsc.pergene.filter01.pi05.n50.bignormal.data.rds
|           |-- humantungipsc.pergene.filter01.pi09.n50.bignormal.data.rds
|       |-- mousekleinesc/
|           |-- mousekleinesc.allgenes.filter01.pi05.n50.bignormal.data.rds
|           |-- mousekleinesc.allgenes.filter01.pi09.n50.bignormal.data.rds
|           |-- mousekleinesc.pergene.filter01.pi05.n50.bignormal.data.rds
|           |-- mousekleinesc.pergene.filter01.pi09.n50.bignormal.data.rds
|       |-- gtex001/
|           |-- gtex001.allgenes.filter01.pi05.n50.bignormal.data.rds
|           |-- gtex001.allgenes.filter01.pi09.n50.bignormal.data.rds
|           |-- gtex001.pergene.filter01.pi05.n50.bignormal.data.rds
|           |-- gtex001.pergene.filter01.pi09.n50.bignormal.data.rds




workflow structures:
|-- submit-makeData.sh
|   |-- makeData.sbatch
|   |-- makeData.R: run ashbun::simulationWrapper()


# extract count table
eset <- get(data("MouseZeiselBrain"))
counts <- exprs(eset)

# filter samples with zero depth
samples_any_gene <- which( colSums(counts) != 0)
counts <- counts[,samples_any_gene]

# filter genes with zero count
genes_to_include <- which(apply(counts,1,sum)>0)
counts <- counts[genes_to_include, ]


# --- register parallel computing clusters
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)

simdata_list <- simulationWrapper(counts, Nsim = 5,
                                  Ngenes = 1000,
                                  Nsamples = 50,
                                  sample_method = "all_genes",
                                  pi0 = .9,
                                  beta_args = args.big_normal(betapi = 1,
                                                              betamu = 0, betasd = .8))

output <- vector("list", 5)
for (index in 1:length(eval_output)) {
  output[[index]]  <- query.evaluation(counts = simdata_list[[index]]$counts,
                                       condition = simdata_list[[index]]$condition,
                                       is_nullgene = simdata_list[[index]]$is_nullgene,
                                       methodsNormalize = c("LIB", "TMM", "census","scran"),
                                       methodsMeanExpression = c("BPSC", "MAST", "ROTS", "SCDE",
                                                                 "DESeq2", "edgeR", "limmaVoom"),
                                       nsim = index)
}

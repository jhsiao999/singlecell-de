library(Biobase)
library(ashbun)
library(singleCellRNASeqMouseZeiselBrain)

# --- Save processed/cleaned raw data
|-- /project2/gilad/joycehsiao/singlecell-de
|    |-- data/
|        |-- mousezeiselbrain.rds
|        |-- humantungipsc.rds
|        |-- mousekleinesc.rds
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


# ---- Make simulated data and save in
|-- /scratch/midway2/joycehsiao/singlecell-de
|   |-- simulated_data/
|       |-- mousezeiselbrain/
|           |-- mousezeiselbrain.allgenes.pi05.n50.bignormal.data.rds
|           |-- mousezeiselbrain.allgenes.pi09.n50.bignormal.data.rds
|           |-- mousezeiselbrain.pergene.pi05.n50.bignormal.data.rds
|           |-- mousezeiselbrain.pergene.pi09.n50.bignormal.data.rds
|       |-- humantungipsc/
|           |-- humantungipsc.allgenes.pi05.n50.bignormal.data.rds
|           |-- humantungipsc.allgenes.pi09.n50.bignormal.data.rds
|           |-- humantungipsc.pergene.pi05.n50.bignormal.data.rds
|           |-- humantungipsc.pergene.pi09.n50.bignormal.data.rds
|       |-- mousekleinesc/
|           |-- mousekleinesc.allgenes.pi05.n50.bignormal.data.rds
|           |-- mousekleinesc.allgenes.pi09.n50.bignormal.data.rds
|           |-- mousekleinesc.pergene.pi05.n50.bignormal.data.rds
|           |-- mousekleinesc.pergene.pi09.n50.bignormal.data.rds



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

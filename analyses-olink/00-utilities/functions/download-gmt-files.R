dir.create("gmt", showWarnings = FALSE)

download.file("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020",
              "../../../data/gmt/MSigDB_Hallmark_2020.gmt", mode = "wb")

download.file("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_Pathways_2024",
              "../../../data/gmt/Reactome_Pathways_2024.gmt", mode = "wb")

download.file("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
              "../../../data/gmt/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt", mode = "wb")

download.file("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human",
              "../../../data/gmt/KEGG_2021_Human.gmt", mode = "wb")
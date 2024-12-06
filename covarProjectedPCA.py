# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 16:25:17 2024

@author: PEIXOTT
"""
import os
import argparse
import numpy as np


def execute(command, logFile=""):
    #os.system(command)
    print(command)



def convertPfileToBfile(pfile, name, folder, plink2, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink2} --pfile {pfile} --make-bed --out {outputPrefix} --double-id --keep-allele-order"

    execute(commandLine)

    return outputPrefix


def mergeCommonFiles(file1, file2, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {file1} --bmerge {file2} --make-bed --out {outputPrefix}"

    execute(commandLine)
    return outputPrefix

def extractVariants(bfile, extractFile, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {bfile} --make-bed --out {outputPrefix} --extract {extractFile}"

    execute(commandLine, logFile)

    return outputPrefix

def removeLDAndMAFToPCA(merged, target, ref, folder, name, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    command = f"{plink1} --bfile {merged} --indep-pairwise 200 50 0.2 --out {folder}/{name} --maf 0.01"
    execute(command, logFile)

    targetLD = extractVariants(target, f"{folder}/{name}.prune.in", f"{name}_Target", folder, plink1, logFile)
    if ref == "":
        return targetLD, ""

    refLD = extractVariants(ref, f"{folder}/{name}.prune.in", f"{name}_Ref", folder, plink1, logFile)

    return targetLD, refLD

def mergeRefAndTarget(bfileTarget, bfileRef, folder, name, plink1, logFile):
    print(f"Merging reference and target")

    execute(f"mkdir {folder}/{name}_MergeRefAlt")

    command = f"{plink1} --bfile {bfileTarget} --make-bed --out {folder}/{name}_MergeRefAlt/TargetMAF --maf 0.01"
    execute(command, logFile)
    execute(f"cp {folder}/{name}_MergeRefAlt/TargetMAF.bed {folder}/{name}_MergeRefAlt/TargetMAFID.bed")
    execute(f"cp {folder}/{name}_MergeRefAlt/TargetMAF.fam {folder}/{name}_MergeRefAlt/TargetMAFID.fam")


    execute(f"cp {bfileRef}.bed {folder}/{name}_MergeRefAlt/Ref.bed")
    execute(f"cp {bfileRef}.fam {folder}/{name}_MergeRefAlt/Ref.fam")


    print(f"Changing the bim ID for target file")
    #Open bim file to change varID to chrom:pos:A1:A2
    dictVar = {}
    bimFileIn = open(f"{folder}/{name}_MergeRefAlt/TargetMAF.bim")
    bimFileOut = open(f"{folder}/{name}_MergeRefAlt/TargetMAFID.bim", "w")
    for line in bimFileIn:
        chrom, ID, poscM, posBP, A1, A2 = line.strip().split()
        newID = f"{chrom}:{posBP}:{A1}:{A2}"
        bimFileOut.write(f"{chrom}\t{newID}\t{poscM}\t{posBP}\t{A1}\t{A2}\n")

        if newID not in dictVar:
            dictVar[newID] = True
        else:
            dictVar[newID] = False
    bimFileIn.close()
    bimFileOut.close()



    print(f"Changing the bim ID for ref file and generating the list of variants in common")
    # Open bim file to change varID to chrom:pos
    bimFileIn = open(f"{bfileRef}.bim")

    bimFileOut = open(f"{folder}/{name}_MergeRefAlt/Ref.bim", "w")
    inCommon = open(f"{folder}/{name}_MergeRefAlt/inCommon.txt", "w")
    for line in bimFileIn:
        chrom, ID, poscM, posBP, A1, A2 = line.strip().split()
        newID = f"{chrom}:{posBP}:{A1}:{A2}"
        bimFileOut.write(f"{chrom}\t{newID}\t{poscM}\t{posBP}\t{A1}\t{A2}\n")

        if newID in dictVar:
            #print(f"In common {newID}")
            if dictVar[newID]:
                #print(f"not added because it was false")
                inCommon.write(f"{newID}\n")

    inCommon.close()
    bimFileIn.close()
    bimFileOut.close()

    refFile = f"{folder}/{name}_MergeRefAlt/Ref"
    targetFile = f"{folder}/{name}_MergeRefAlt/TargetMAFID"
    commonList = f"{folder}/{name}_MergeRefAlt/inCommon.txt"
    targetCommon = extractVariants(targetFile, commonList, f"Target_common", f"{folder}/{name}_MergeRefAlt/", plink1, logFile)
    refCommon = extractVariants(refFile, commonList, f"Ref_common", f"{folder}/{name}_MergeRefAlt/", plink1, logFile)

    return mergeCommonFiles(targetCommon, refCommon, "Merged", f"{folder}/{name}_MergeRefAlt/", plink1, logFile), targetCommon, refCommon


def getProjectedPCA(target, reference, gcta, X, threads, name, folder, plink1, logFile):
    command = f"{gcta} --bfile {reference} --maf 0.01 --make-grm --out {folder}/{name} --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile)

    command = f"{gcta} --grm {folder}/{name} --maf 0.01 --pca 50 --out {folder}/{name}_PCs --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile)

    command = f"{gcta} --bfile {reference} --maf 0.01 --pc-loading {folder}/{name}_PCs --out {folder}/{name}_VariantLoading --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile)
    if target == "":
        return f"{folder}/{name}_PCs"


    command = f"{gcta} --bfile {target} --maf 0.01 --project-loading {folder}/{name}_VariantLoading 50 --out {folder}/{name}_TargetPC --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile)

    return f"{folder}/{name}_TargetPC"


#======================================================================== Covar =========================================================================
def modelSelection(covar, modelSelection, folder, name, logFile):
    commandLine = f"Rscript {modelSelection} {covar} {folder}/{name}"

    execute(commandLine)

    print(f"Opening the file {folder}/{name}_variables.tsv")
    fileModel = open(f"{folder}/{name}_variables.tsv")
    model = []
    for line in fileModel:
        model.append(line.strip())

    print(f"Return the model: {model}")
    return model



def buildCovarFile(covarDict, outlier, PCList, PCASource, bfile, folder, name):
    outlierList = []

    #outlierFile = open(outlier)
    #for line in outlierFile:
    #    IID, FID = line.strip().split()
    #    outlierList.append(IID)
    #outlierFile.close()

    covarList = []
    famFile = open(f"{bfile}.fam")
    covarFile = open(f"{folder}/{name}.tsv", "w")

    #Header -> IID <covar non PC> <Covar PC in Source_PCNumber format>
    covarFile.write("IID")
    for ind in covarDict:
        for covar in covarDict[ind]:
            if covar != "PCA":
                covarFile.write(f"\t{covar}")
                covarList.append(covar)
        break
    for source in PCASource:
        for PC in PCList:
            covarFile.write(f"\t{source}_{PC}")
    covarFile.write("\n")

    for line in famFile:
        FID, IID, mother, father, sex, pheno = line.strip().split()
        if IID not in outlierList and IID in covarDict:
            covarFile.write(f"{IID}")
            for covar in covarList:
                covarFile.write(f"\t{covarDict[IID][covar]}")
            for source in PCASource:
                for PC in PCList:
                    covarFile.write(f"\t{covarDict[IID]['PCA'][source][PC]}")
            covarFile.write("\n")
    covarFile.close()
    famFile.close()
    return f"{folder}/{name}.tsv"


def createOutlierList(covarDict, bfile, PCASource, PCList, folder, name, missing, model):
    print(f"{bfile}.fam")
    #print(f"{covarDict}")
    famFile = open(f"{bfile}.fam")
    sampleList = []

    for line in famFile:
        FID, IID, mother, father, sex, pheno = line.strip().split()
        sampleList.append(IID)

    removalDict = {}

    for dataSource in PCASource:
        
        #Get all PCs that are in the model
        for PC in PCList:
            if PC in model:
                listPC = []
                
                #Get all the samples that are on covarList
                for ind in sampleList:
                    if ind in covarDict:
                        PCInfo = float(covarDict[ind]["PCA"][dataSource][PC])
                        listPC.append(PCInfo)
    
                mean = np.mean(listPC)
                sd = np.std(listPC)
    
                lowerBound = mean - 3 * sd
                upperBound = mean + 3 * sd
    
    
                for ind in sampleList:
                    if ind in covarDict:
                        PCInfo = float(covarDict[ind]["PCA"][dataSource][PC])
                        if PCInfo < lowerBound or PCInfo > upperBound:
                            if ind not in removalDict:
                                removalDict[ind] = []
                                
                            print(f"{ind} outlier {PC}")
                            removalDict[ind].append(f"{PC}({dataSource})")

    fileToRemove = open(f"{folder}/{name}_outlierList.txt", "w")
    fileToRemoveExplanation = open(f"{folder}/{name}_outlierExplanation.txt", "w")



    for ind in missing:
        fileToRemove.write(f"{ind}\t{ind}\n")
        fileToRemoveExplanation.write(f"{ind}: missing data\n")
    
    for ind in removalDict:
        fileToRemove.write(f"{ind}\t{ind}\n")
        fileToRemoveExplanation.write(f"{ind}:")
        for info in removalDict[ind]:
            fileToRemoveExplanation.write(f" {info}")
        fileToRemoveExplanation.write(f"\n")
        missing.append(ind)
    
    fileToRemove.close()
    fileToRemoveExplanation.close()

    return f"{folder}/{name}_outlierList.txt", missing



def addPCAToCovarDict(filePCA, covarDict, dataSource):
    for ind in covarDict:
        #print(f"For error: {covarDict[ind]}")
        if "PCA" not in covarDict[ind]:
            covarDict[ind]["PCA"] = {}
        if dataSource not in covarDict[ind]["PCA"]:
            covarDict[ind]["PCA"][dataSource] = {}


    file = open(f"{filePCA}.proj.eigenvec")

    #Project PCA has no header, while the previous PCA had
    PCList = []

    for line in file:
        data = line.strip().split("\t")
        ind = data[1]
        if ind in covarDict:
            for i in range(2, len(data)):
                PC = f"PC{i-1}"
                if PC not in PCList:
                    PCList.append(PC)
                covarDict[ind]["PCA"][dataSource][PC] = data[i]


    return covarDict, PCList


def readCovarFile(covarTable):
    print('We are reading the covar table. We are assuming that the ID is the first col')
    print('We are also assuming that there is the column SEX and the Phenotype column is named DISEASE')
    print('We are also checking if there is any covar field that is empty or NA')

    file = open(covarTable)

    covarDict = {}
    missing = []
    header = True

    # doNotFix is a flag to put the phenotype in PLINK2 format (1 control, 2 case)
    doNotFix = False
    for line in file:
        if header:
            header = False
            splitHeader = line.strip().split()
        else:
            split = line.strip().split()

            nonNA = True
            for i in range(len(split)):
                if split[i] == "NA" or split[i] == "" or split[i] == " " or split[i] == "nan":
                    nonNA = False
                    ind = split[0]
                    data = split[i]
                    headerName = splitHeader[i]

                    print(f'Removing the ind {ind} because there is missing data ({data}) on the field {headerName}')
                    missing.append(ind)

            if nonNA:
                covarDict[split[0]] = {}
                for i in range(1, len(split)):
                    if splitHeader[i].upper() == "DISEASE":
                        if split[i] == "0" or split[i] == 0:
                            split[i] = "1"
                        elif split[i] == "1" or split[i] == 1:
                            split[i] = "2"
                        elif split[i] == "2" or split[i] == 2:
                            split[i] = "3"
                            doNotFix = True
                        else:
                            print(f"Unknown status value to sample {split[0]}: {split[i]} ")
                    covarDict[split[0]][splitHeader[i].upper()] = split[i]

    if doNotFix:
        for sample in covarDict:
            if covarDict[sample]["DISEASE"] == "2":
                covarDict[sample]["DISEASE"] = "1"
            elif covarDict[sample]["DISEASE"] == "3":
                covarDict[sample]["DISEASE"] = "2"
            else:
                print(f"Unknown sample DISEASE field:  {covarDict[sample]['DISEASE']}")

    return covarDict, missing

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PCA and regression')

    data = parser.add_argument_group("Data arguments")
    data.add_argument('-A', '--autosomal', help='Genotyped file name for autosomal chromosome', required=True)
    data.add_argument('-t', '--tableCovar', help='File with covariatives to be added to the model', required=True)
    data.add_argument('--threads', help='Number of processors to be used (default = 1)',
                      required=False, default = 1)

    refData = parser.add_argument_group("Reference data arguments")
    refData.add_argument('-R', '--AutosomalRef', help='Data with parental reference data to run autosomal PCA (optional)', required=False, default="")

    output = parser.add_argument_group("Output arguments")
    output.add_argument('-n', '--name', help='Analysis name', required=True)
    output.add_argument('-f', '--folder', help='Folder to output files', required=True)

    pcaArg = parser.add_argument_group("PCA arguments")
    pcaArg.add_argument('--model', help='Regression model. If you do not provide a model the script will'
                                        'use the selectModel.R to build the model', required=False, default="", nargs="+")


    programs = parser.add_argument_group("Programs")
    programs.add_argument('--plink2', help='Path of PLINK 2 (default = plink2)', required=False, default="plink2")
    programs.add_argument('--plink1', help='Path of PLINK 1 (default = plink)', required=False, default="plink")
    programs.add_argument('--gcta', help='Path of gcta', required=True)
    programs.add_argument('--selectModel', help='Path of selectModel script', required=False)

    args = parser.parse_args()
    logFile = ""

    covarDict, missing = readCovarFile(args.tableCovar)
    os.system(f"mkdir {args.folder}")
    
    if os.path.isfile(f"{args.autosomal}.pgen"):
        args.autosomal = convertPfileToBfile(args.autosomal, "Target", args.folder, args.plink2, logFile)


    if args.AutosomalRef != "":
        # Projected PCA
        if os.path.isfile(f"{args.AutosomalRef}.pgen"):
            args.AutosomalRef = convertPfileToBfile(args.AutosomalRef, "Ref", args.folder, args.plink2, logFile)


        bfileMerged, targetCommon, refCommon = mergeRefAndTarget(args.autosomal, args.AutosomalRef, args.folder, f"{args.name}_AutosomalPCA", args.plink1, logFile)
        targetLD, refLD = removeLDAndMAFToPCA(bfileMerged, targetCommon, refCommon, args.folder, f"{args.name}_LD", args.plink1, logFile)
        autosomalPCA = getProjectedPCA(targetLD, refLD, args.gcta, False, args.threads, "AutosomalPCA", args.folder, args.plink1, logFile)
    else:
        #PCA without referece
        #Send ref as empty to cancel the LD removal for ref
        targetLD, refLD = removeLDAndMAFToPCA(args.autosomal , args.autosomal, "", args.folder, f"{args.name}_LD",args.plink1, logFile)

        #Send target as empty to perform PCA only on target (that will be the ref on the function)
        autosomalPCA = getProjectedPCA("", targetLD, args.gcta, False, args.threads, "AutosomalPCA",
                                       args.folder, args.plink1, logFile)

    covarDict, PCList = addPCAToCovarDict(autosomalPCA, covarDict, "GCTA")
    
    PCASource = ["GCTA"]
    outlierList = []
    
    covar = buildCovarFile(covarDict, outlierList, PCList, PCASource, args.autosomal, args.folder, f"{args.name}")
    model = modelSelection(covar, args.selectModel, args.folder, f"{args.name}", logFile)
    outlierListFile, outlierList = createOutlierList(covarDict, args.autosomal, PCASource, PCList, args.folder, f"{args.name}", missing, model)

    covarAll = []
    covarOrder = []
    
    print(covarDict)
    
    for ind in covarDict:
        for covar in covarDict[ind]:
            if covar == "PCA":
                for PC in covarDict[ind][covar]["GCTA"]:
                    PCName = f"GCTA_{PC}"
                    if PCName not in covarAll:
                        covarAll.append(PCName)
            else:
                if covar not in covarAll:
                    covarAll.append(covar)
    print(f"All covar: {covarAll}")

    

    fileFinal = open("CovarPCAProjected.tsv", "w")
    covarOrder.append("DISEASE")
    fileFinal.write(f"IID\tDISEASE")

    for covar in covarAll:
        if covar in model:
            fileFinal.write(f"\t{covar}")
            print(f"The covar {covar} is on the model")
            covarOrder.append(covar)
        else:
            for modelCovar in model:
                if covar in modelCovar and "PC" not in covar:
                    fileFinal.write(f"\t{covar}")
                    print(f"The covar {covar} is on the model ({modelCovar})")
                    covarOrder.append(covar)
                    break
  
    fileFinal.write(f"\n")
    for ind in covarDict: 
        print(f"Addind {ind} -> ", end = "")
        if ind not in outlierList:
            print(f"Addind")
            if "PC1" in covarDict[ind]['PCA']['GCTA']:
                fileFinal.write(f"{ind}")
                for covar in covarOrder:
                    print(f"{covar}")    
                    if "PC" in covar:
                        PC = covar.split("_")[-1]
                        fileFinal.write(f"\t{covarDict[ind]['PCA']['GCTA'][PC]}")
                    else:
                        fileFinal.write(f"\t{covarDict[ind][covar]}")
                fileFinal.write(f"\n")
        else:
            print(f"Not in")
    fileFinal.close()

%% Epigenetic Analysis Pipeline

% This repository contains the analysis pipeline and scripts for the study
% "Mapping Epigenetic Gene Variant Dynamics: Comparative Analysis of
% Frequency, Functional Impact, and Trait Associations in African and
% European Populations." It includes custom code for data processing,
% statistical analysis, and visualization, aiming to provide insights into
% the genetic and epigenetic differences between African and European
% populations. The processed data will be available in a Zenodo repository.

clc; clear; close all force

% If you working on a computer with less than 60GB on ram you may need to
% set the option below to false
multiCoreHighMem = false; % true -- for complete excution of pipeline

% add the path the contains the files
addpath(genpath('/required_files'))

%% First get the epigenetic genes

% here are the epigenetic genes taken from reactome pathway database and
% the supplementary file of: Gnad et al. BMC Genomics 2015, 16(Suppl 8):S5
% http://www.biomedcentral.com/1471-2164/16/S8/S5. Check the
% epigeneticPipeline.m for more details
epigeneticGenes = readtable('epigeneticsGenes.xlsx');
epigeneticGenes.Properties.VariableNames = ["HugoSymbol","GeneClass"];

% get the counts of the epigenetic genes and sort the table
epigeneticGenes.GeneClass = categorical( epigeneticGenes.GeneClass) ;
geneCounts = table(categories(epigeneticGenes.GeneClass), ...
    countcats(epigeneticGenes.GeneClass), 'VariableNames', ...
    {'GeneClass','Count'}) ;

% arrange by descending order and convert to categorical
geneCounts = sortrows(geneCounts,'Count','ascend') ;
geneCounts.GeneClass = categorical(geneCounts.GeneClass,...
    geneCounts.GeneClass) ;

% save the genes to the supplementary data 
writetable( epigeneticGenes, 'Supplementary Data 1.xlsx', ...
    'Sheet','Epigenetic Gene')

% here are the high_quality variants
high_qual_variants = readtable( ...
    'epigenetic_gene_high_quality_variants.txt');
high_qual_variants = high_qual_variants( ...
    high_qual_variants.af_AFR > 1e-5 &  ...
    high_qual_variants.af_EUR > 1e-5 , :) ;

% remove the palette variable
clear palette

%% Distribution of SNPs

% What is the distribution of epigenetic-related SNPs across different
% populations - African and Europeans and across different African regions?
% - UK Biobank data (Europeans). 

if ~exist('imputed_snp_freq.csv','file')
    
    fprintf('\Now getting the epigenetic gene frequencies\n')
    
    % I have already processesed the SNPs on the cluster and here is the
    % file
    snpFreq = readtable('UKB_SNP_Comparisons_Less.csv') ;
    
    % return only the epigenetic genes frequencies
    snpFreq = innerjoin(epigeneticGenes,snpFreq) ;
    
    % save to a CSV file
    writetable(snpFreq,'imputed_snp_freq.csv')
else
    % load the processed datasets of only snpfreqs of epigenetic genes
    snpFreq = readtable('imputed_snp_freq.csv') ;
    
    % add the sig groups for each genes in the snp freq data
    snpFreq.SigCount = snpFreq.FDR < 0.05 ;
end

% get the variant matrix for only epigenetic genes
if ~exist('variantqc_epi.txt','file')

    fprintf('\nCreating a new variantqc for epigenetic\n')
    % load the variant qc data
    variantqc = readtable('variant_qc_metrics.txt') ;

    % merge the data with epigenetic genes data
    variant_epi = innerjoin(variantqc, ...
        snpFreq(:,{'Variant','HugoSymbol'}), ...
        "LeftKey","rsid","RightKey","Variant") ;

    % save the file
    writetable(variant_epi,'variantqc_epi.txt')

    % exit the code 
    return
end

% ***********************************************************************
% eQTLs, hQTLs, mQTL and sQTLs   
% ***********************************************************************

% get the mQTLs from the cluster 
if ~exist('epigenetic_enhancer.txt','file')
    addpath('/required_files/oncoBase')

    fprintf('\n Processing enhancer data\n')
    % laod the enhance data and return only those of epigenetic genes
    enhancer = readtable('JEME.txt') ;
    enhancer.Properties.VariableNames(4) = "HugoSymbol" ;
    theCommon = ismember(enhancer.HugoSymbol, snpFreq.HugoSymbol);
    enhancer = enhancer(theCommon,:) ;
    writetable(enhancer,'epigenetic_enhancer.txt') ;

    fprintf('\n Processing mQTL data\n')
    % load the mQTL data and return only those for epigenetic genes
    mqtl = readtable('mQTL.txt');
    mqtl.Properties.VariableNames(2) = "Variant" ;
    mqtl = innerjoin(snpFreq, mqtl ,"Keys","Variant")  ;
    writetable(mqtl,'epigenetic_mQTL.txt') ;

    fprintf('\n Processing hQTL data\n')
    % load the hQTL data and return only those for epigenetic genes
    % histone quantitative trait loci (hQTLs)
    hqtl = readtable('local_hQTLs.tab.txt');
    hqtl.Properties.VariableNames(2) = "Variant" ;
    hqtl = innerjoin(snpFreq, hqtl ,"Keys","Variant")  ;
    writetable(hqtl,'epigenetic_hQTL.txt') ;

    fprintf('\nThe mQTL, hQTL and enhancer data has been procesed\n')

    % end the program 
    % return 

    % change back the direct
    cd('..')
end

% process the sqtls data
if ~exist('epigenetic_sQTL.txt','file')

    % set up an empty matrix
    sqtl = [] ;

    % define the tissue types available that are aviable in the folder
    if exist('/required_files/oncoBase/sQTLs.GTEx.V7.RSEM','dir')
        addpath('/required_files/oncoBase/sQTLs.GTEx.V7.RSEM')
        cd('/required_files/oncoBase/sQTLs.GTEx.V7.RSEM')
        listing = ...
            dir('/required_files/oncoBase/sQTLs.GTEx.V7.RSEM');
        listing = struct2table(listing) ;
        listing = listing.name( ~ismember(listing.name,{'.','..'})) ;
    end

    for ii = 1:length(listing)
        % print something to the screen
        fprintf('\nGetting sQTLs of %s \n', listing{ii})

        % change the directory to the tissue folder
        if isfolder(listing{ii})
            cd(listing{ii})
        end

        % unzip the file the file if this is not already done
        if exist('sqtls-0.05fdr.permuted.tsv.gz','file')
            gunzip('sqtls-0.05fdr.permuted.tsv.gz')
        end

        % get the current sQTL data and return only the requried columns
        cur_sqtl = readtable('sqtls-0.05fdr.permuted.tsv', ...
            'FileType','text');
        cur_sqtl = cur_sqtl(:,{'geneId','snpId','fdr'}) ;
        cur_sqtl = unique(cur_sqtl) ;
        cur_sqtl = cur_sqtl(cur_sqtl.fdr < 0.05 ,:) ;
        cur_sqtl = addvars(cur_sqtl, ...
            repmat( listing(ii), height(cur_sqtl), 1) , ...
            'NewVariableNames',{'tissue'},'Before',1) ;

        % get the require columns of the data
        sqtl = [sqtl; cur_sqtl] ;
    end
    
    % the gene ids in the eQTLs are ENSG gene ids so I have to convert them
    % to hugo gene symbols
    try
        allBiomart = readtable('mart_export.txt') ;
    catch
        allBiomart = readtable('/required_files/mart_export.txt');
    end

    allBiomart = allBiomart( ...
        ismember(allBiomart.HGNCSymbol,epigeneticGenes.HugoSymbol),: ) ;
    allBiomart = allBiomart(:,{'GeneStableIDVersion','HGNCSymbol'}) ;
    allBiomart.Properties.VariableNames = {'geneId','HugoSymbol'} ;
    sqtl = innerjoin(sqtl, allBiomart);

    % clean up the data
    sqtl.Chrom = extractBefore(sqtl.snpId,'_') ;
    sqtl.Positions = extractBetween(sqtl.snpId,'_','_') ;
    sqtl.geneId = [];
    sqtl = splitvars(sqtl) ;
    sqtl.Positions_2 = [];
    sqtl.Properties.VariableNames(end) = "Position" ;

    % make the data workable with future fuctions
    sqtl.Position = str2double(sqtl.Position);
    sqtl.Chrom = str2double(sqtl.Chrom) ;

    % change back the direct
    cd('..')

    % save the file
    writetable(sqtl,'epigenetic_sQTL.txt') ;

    clear allBiomart cur_sqtl ii listing
end

% convert the genes and sig group to categorical
snpFreq.HugoSymbol = categorical(snpFreq.HugoSymbol) ;
snpFreq.highGroup = categorical(snpFreq.highGroup);
snpFreq.GeneClass = categorical(snpFreq.GeneClass);

% ************** Return only the High Confidence SNPs ***************
fprintf('\nThe # of all variants is %d\n', height(snpFreq) )
snpFreq = snpFreq( ismember(snpFreq.Variant, high_qual_variants.rsid), :);
fprintf('\nThe # of high-quality high-freq variants is %d\n', ...
    height(snpFreq) )
% *******************************************************************

% add information of the gene length to the snp freq table using
% infomration from biomart
biomart = readtable('mart_export_with_TSS_ALL_GENES.txt') ;

% get the length of the genes
biomart.GeneLength = biomart.GeneEnd_bp_ - biomart.GeneStart_bp_ ;
biomart = biomart(:,{'GeneName','GeneLength'});

% get the unique genes and convert the gene name to hugo symbol
[~,uniqueRows]= unique(biomart.GeneName) ;
biomart = biomart(uniqueRows,:);
biomart.Properties.VariableNames(1) = "HugoSymbol" ;

% add the epigenetic genes class to the biomart data
biomart = innerjoin(biomart,epigeneticGenes) ;

% convert to Hugo Symbol
biomart.HugoSymbol = categorical(biomart.HugoSymbol) ;

% here is the summary of snps that vary significantly between groups
summary(snpFreq.highGroup)

% here is the number of samples
ukbAFRsamples = 5978  ;
ukbEURsamples = 383471;

% What epigenetic genes and gene classes has the most variants and also the
% most significantly different variants between Africans and Europeans

% here are the vars for the bar plots
theVars = ["GeneClass","HugoSymbol"] ;
theColors = [0 0.451 0.761; 0.761 0 0.451] ;

% *********** plot a the frist bar graph of epigenetic genes **********
figure()
set(gcf,'position',[500,500,1200,700]);
tiledlayout(2,3,'TileSpacing','compact')
clf ;

% here the letters for the plots
theLetters = 'a':'j';

% ***************** The bar graph of epigenetic genes ***************
nexttile

barh(geneCounts.GeneClass,geneCounts.Count)

% edit the axis and and adjust the figure
xlabel('Number of genes')
set(gca,'FontSize',12,'LineWidth',1,'Box','off','TickDir','out')

% add a letter to the figure and remove the letter from the array
text(-0.2, 1 ,theLetters(1),'Units','normalized','FontWeight','bold', ...
    'FontSize',28)
theLetters(1) = [];

% *********************************************************************

% ******************* the scatter plot of epigenetic genes **********
% populations Take all populations and group them into African and
% non-African -select cut-off frequency as definition for present/absent
% (remove admixed pops) Share with Akin and Gaone: SNP frequencies of
% Africans and Europeans (UK Biobank). -

% the colors to use the for two groups
groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10  ] ;

% plot the graph before removing the causal SNPs
nexttile
hold on

% plot for the EUR
scatter(snpFreq.af_AFR(snpFreq.highGroup == "EUR"), ...
    snpFreq.af_EUR(snpFreq.highGroup == "EUR"), 7,'filled', ...
    'MarkerFaceColor',groupColors(2,:),'Marker','o',...
    'MarkerEdgeColor',groupColors(2,:),'MarkerFaceAlpha',0.7)

% for the AFR
scatter(snpFreq.af_AFR(snpFreq.highGroup == "AFR"), ...
    snpFreq.af_EUR(snpFreq.highGroup == "AFR"), 7,'filled', ...
    'MarkerFaceColor',groupColors(1,:),'Marker','o',...
    'MarkerEdgeColor',groupColors(1,:))

% plot the non significant data
nonSigFisher = snpFreq.FDR > 0.05 ;
scatter(snpFreq.af_AFR(nonSigFisher), ...
    snpFreq.af_EUR(nonSigFisher), ...
    7,'filled','MarkerFaceColor',[ 0.5 0.5 0.5 ],'Marker','o',...
    'MarkerEdgeColor',[ 0.5 0.5 0.5 ],'MarkerFaceAlpha',0.7)

% edit the chart elements
set(gca,'FontSize',12,'LineWidth',1,'Box','off','TickDir','out')
ylabel('SNP Frequency in EUR')
xlabel('SNP Frequency in AFR');
legend({sprintf('EUR: %d SNPs',sum(snpFreq.highGroup=="EUR")),...
    sprintf('AFR: %d SNPs',sum(snpFreq.highGroup=="AFR")),...
    sprintf('Non Sig: %d SNPs',sum(snpFreq.highGroup=="None")),...
    },'Location','best')
% title('Frequency Difference','FontSize',16,'FontWeight','bold')

% add a letter to the figure and remove the letter from the array
text(-0.2, 1 ,theLetters(1),'Units','normalized','FontWeight','bold', ...
    'FontSize',28)
theLetters(1) = [];

% *******************************************************************

% plot the figures in a loop
for ii = 1:length(theVars)

    % the variable for normalising the SNPs count with the length of the
    % genes

    % get the data for barplots which has the snpfreqs for all genes
    curData = groupsummary(snpFreq, theVars{ii}, "sum","SigCount") ;

    % sort the table
    if ii == 1
        curData = sortrows(curData,'GroupCount','ascend') ;
    else
        curData = sortrows(curData,'GroupCount','descend') ;
    end

    % get the top 50 genes
    if height(curData) > 50
        curData = curData(1:50,:);
    end

    % rearrage the variables
    curData.(theVars{ii}) = categorical( curData.(theVars{ii}), ...
        curData.(theVars{ii}) );
    
    fprintf('\nHere is the data for %s\n',theVars{ii})
    disp(curData)

    % create a new figure only the second plots
    if strcmp(theVars{ii},"GeneClass" )

        % here is the next figure
        nexttile
        % here is the bargraph for the mean values
        barh(curData.(theVars{ii}) ,curData.sum_SigCount,...
            'FaceColor',theColors(ii,:), 'EdgeColor', theColors(ii,:));

        % add for the max values
        hold on
        barh(curData.(theVars{ii}),curData.GroupCount, ...
            'FaceColor', theColors(ii,:), 'EdgeColor', ...
            theColors(ii,:),'FaceAlpha',0.3,'EdgeAlpha',0.3);

        % edit the axis and and adjust the figure
        xlabel('Number of SNPs')
        set(gca,'FontSize',12,'LineWidth',1,'Box','off','TickDir','out')

        % convert the axis labels
        b = get(gca,'XTick')/1e3 ;
        b = split( cellstr(num2str(b)) ) ;
        b =  strcat(b, 'K')';
        set(gca,'XTickLabel',b )

        % add a letter to the figure and remove the letter from the array
        text(-0.2, 1 ,theLetters(1),'Units','normalized', ...
            'FontWeight','bold', 'FontSize',28)
        theLetters(1) = [];

    else

        nexttile([1,3])
        % here is the bargraph for the mean values
        bar(curData.(theVars{ii}),curData.sum_SigCount, ...
            'FaceColor',theColors(ii,:),'EdgeColor', theColors(ii,:));

        % add for the max values
        hold on
        bar(curData.(theVars{ii}),curData.GroupCount, ...
            'FaceColor', theColors(ii,:),'EdgeColor', theColors(ii,:),...
            'FaceAlpha',0.3,'EdgeAlpha',0.3);

        % edit the axis and and adjust the figure
        ylabel('Number of SNPs')
        set(gca,'FontSize',12,'LineWidth',1,'Box','off',...
            'XTickLabelRotation',90,'TickDir','out')

        % put a normalise on the y-axis
        % convert the axis labels
        b = get(gca,'YTick')/1e3 ;
        b = strcat( cellstr(num2str(b')) ,'K');
        set(gca,'YTickLabel',b )

        % change the title of the figures based on the type of data
        xlabel('Genes') ;

        % add a letter to the figure and remove the letter from the array
        text(-0.05, 1 ,theLetters(1),'Units','normalized', ...
            'FontWeight','bold', 'FontSize',28)
        theLetters(1) = [];
    end

    % add a legend and title to the figure
    legend({'Sig SNPs','Total SNPs'},'Box','off','Location','best')
end

hold off

clear ii curData jj curBiomart uniqueRows b

% save a clean version of snpFreq to used in teh
% final_ukb_locusZoomPlot_trait.m file 
writetable(snpFreq,'highConf_snpFreq.txt')

%% Plot the figures for the normalised genes

% *********** plot a the frist bar graph of epigenetic genes **********
tiledlayout(1,3)

% here the letters for the plots
theLetters = 'a':'j';

% plot the figures in a loop
for ii = 1:length(theVars)

    % get the data for barplots which has the snpfreqs for all genes
    curData = groupsummary(snpFreq, theVars{ii}, "sum","SigCount") ;

    % normalise for the length of the coding sequence

    % get the avarage length of the genes class
    if strcmp(theVars{ii},"GeneClass" )
        curBiomart = removevars(biomart,'HugoSymbol') ;

        % get the mean length of gene in each class and remove the
        % unwanted variables
        curBiomart = groupsummary(curBiomart,"GeneClass", ...
            "sum","GeneLength") ;
        curBiomart.Properties.VariableNames(end) = "GeneLength" ;
        curBiomart.GroupCount = [];

    elseif strcmp(theVars{ii},"HugoSymbol" )
        curBiomart = removevars(biomart,'GeneClass') ;
    end

    % Not all the epigenetic genes in the data so only do this for
    % the genes are present
    curData = innerjoin(curData,curBiomart,'Key',theVars{ii}) ;

    % here is the number of samples
    numSamples = sum(ukbAFRsamples,ukbEURsamples) ;

    % divide by the length of the genes
    curData.GroupCount = ...
        (curData.GroupCount./curData.GeneLength)/numSamples;
    curData.sum_SigCount = ...
        ( curData.sum_SigCount./curData.GeneLength )/numSamples;

    % get only the unque rows
    curData = unique(curData) ;

    % sort the table
    curData = sortrows(curData,'GroupCount','descend') ;

    % get the top 50 genes
    if height(curData) > 50
        curData = curData(1:50,:);
    end

    % rearrage the variables
    curData.(theVars{ii}) = categorical( curData.(theVars{ii}), ...
        curData.(theVars{ii}) );

    % here is the figure
    if ii == 1
        nexttile
    else
        nexttile([1,2])
    end

    % for the second figure except for the unnormalised gene class data
    % here is the bargraph for the mean values
    bar(curData.(theVars{ii}),curData.sum_SigCount,'FaceColor', ...
        theColors(ii,:), 'EdgeColor', theColors(ii,:));

    % add for the max values
    hold on
    bar(curData.(theVars{ii}),curData.GroupCount, ...
        'FaceColor', theColors(ii,:), ...
        'EdgeColor', theColors(ii,:),...
        'FaceAlpha',0.3,'EdgeAlpha',0.3);

    % edit the axis and and adjust the figure
    ylabel('Number of SNPs')
    set(gca,'FontSize',12,'LineWidth',1,'Box','off',...
        'XTickLabelRotation',90)

    % put a normalise on the y-axis
    ylabel('Normalised number of SNPs')

    % change the title of the figures based on the type of data
    if strcmp(theVars{ii},"HugoSymbol" )
        % for the normalised snps
        xlabel('Genes') ;
        title('SNPs in Top 50 Genes (normalised)','FontSize',16)
    else
        % for the normalised snps
        xlabel('Gene Classes');
        title('SNPs (normalised)','FontSize',16)
    end

    % add a legend and title to the figure
    legend({'Sig SNPs','Total SNPs'},'Box','off')

    if ii == 1
        % add a letter to the figure and remove the letter from the array
        text(-0.2, 1 ,theLetters(1),'Units','normalized', ...
            'FontWeight','bold', 'FontSize',28)
    else
        % add a letter to the figure and remove the letter from the array
        text(-0.07, 1 ,theLetters(1),'Units','normalized', ...
            'FontWeight','bold', 'FontSize',28)
    end

    theLetters(1) = [];
end

hold off
%% Which SNPs are unique to Africa or skewed in their frequencies 

% populations Take all populations and group them into African and
% non-African -select cut-off frequency as definition for present/absent

% the colors to use the for two groups
groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10  ] ;

% plot the graph before removing the causal SNPs
figure()
hold on

% plot for the EUR
scatter(snpFreq.af_AFR(snpFreq.highGroup == "EUR"), ...
    snpFreq.af_EUR(snpFreq.highGroup == "EUR"), 7,'filled', ...
    'MarkerFaceColor',groupColors(2,:),'Marker','o',...
    'MarkerEdgeColor',groupColors(2,:),'MarkerFaceAlpha',0.7)

% for the AFR
scatter(snpFreq.af_AFR(snpFreq.highGroup == "AFR"), ...
    snpFreq.af_EUR(snpFreq.highGroup == "AFR"), 7,'filled', ...
    'MarkerFaceColor',groupColors(1,:),'Marker','o',...
    'MarkerEdgeColor',groupColors(1,:))

% plot the non significant data
nonSigFisher = snpFreq.FDR > 0.05 ;
scatter(snpFreq.af_AFR(nonSigFisher), ...
    snpFreq.af_EUR(nonSigFisher), ...
    7,'filled','MarkerFaceColor',[ 0.5 0.5 0.5 ],'Marker','o',...
    'MarkerEdgeColor',[ 0.5 0.5 0.5 ],'MarkerFaceAlpha',0.7)

% edit the chart elements
set(gca,'FontSize',14,'LineWidth',1.5,'Box','off')
ylabel('SNP Frequency in EUR')
xlabel('SNP Frequency in AFR');
legend({sprintf('EUR: %d SNPs',sum(snpFreq.highGroup=="EUR")),...
    sprintf('AFR: %d SNPs',sum(snpFreq.highGroup=="AFR")),...
    sprintf('Non Sig: %d SNPs',sum(snpFreq.highGroup=="None")),...
    },'Location','best')
title('Frequency Difference','FontSize',16,'FontWeight','bold')

hold off

% save the figure
saveas(gca,'Frequency_Difference_Scatter.fig','fig')
saveas(gca,'Frequency_Difference_Scatter.png','png')

clear nonSigFisher

%% Plot the Bin Histograp

% Calculate the log of frequencies and handle zeros or negative values if
% necessary.
afr_freqs = snpFreq.af_AFR;
eur_freqs = snpFreq.af_EUR;

% Replace non-positive values with NaN to avoid complex numbers
afr_freqs(afr_freqs <= 0) = nan; 

% Create a binned scatter plot of the data points.
figure();
binscatter(eur_freqs, afr_freqs,'ShowEmptyBins', 'off');
% Set color map to 'jet' for the current axes
colormap(gca, 'parula'); 

% Set color axis scaling based on data range
% caxis([min(afr_freqs, [], 'omitnan'), max(afr_freqs, [], 'omitnan')]); 
colorbar; % Show the color bar for reference

% Edit plot properties.
xlabel("Patient's Age"); 
ylabel('Log Mutation Count');
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'off', ...
    'FontWeight', 'bold'); % , 'YLim', [-1 11]

% Calculate Pearson's linear correlation for complete rows.
[r, pValue] = corr(eur_freqs, afr_freqs, 'Type', 'Pearson', 'Rows', 'complete');

% Annotate the plot with the correlation coefficient and p-value.
% Note: `convertPValue2SuperScript` is a custom function not defined in this code.
% Replace it with the actual function or code to convert p-value to superscript.
annotationText = sprintf("R = %.2f,  P = %.2e", r, pValue);
annotation('textbox', [0.5, 0.95, 0, 0], 'String', annotationText,...
    'FitBoxToText', 'on','FontWeight', 'bold', 'FontSize', 12,...
    'HorizontalAlignment', 'center');

% Add a linear fit line to the plot.
% Handle missing data for the fit.
afr_freqs = fillmissing(afr_freqs, 'nearest');
eur_freqs = fillmissing(eur_freqs, 'nearest');


% Fit and plot.
% Linear regression
b = regress(afr_freqs, [ones(length(eur_freqs), 1), eur_freqs]);
yCalc = [ones(length(eur_freqs), 1), eur_freqs] * b;
hold on; % Keep the existing plot
plot(eur_freqs, yCalc, 'k-', 'LineWidth', 2); % Plot the linear fit
hold off;

%% a scatter hist 

% Assume snpFreq is a table or struct with af_AFR and af_EUR as fields
% for African and European allele frequencies, respectively.

% Calculate the log of frequencies and handle zeros or negative values if
% necessary.
afr_freqs = snpFreq.af_AFR;
eur_freqs = snpFreq.af_EUR;

% Replace non-positive values with NaN to avoid complex numbers
afr_freqs(afr_freqs <= 0) = nan;
eur_freqs(eur_freqs <= 0) = nan;

% % Take the natural log of positive frequencies
% afr_freqs = log(afr_freqs);
% eur_freqs = log(eur_freqs);

% Create a scatter plot with histograms on both axes to visualize the
% relationship between allele frequencies in the AFR and EUR populations.
% Customize the appearance of the plots.
scatterhist(eur_freqs, afr_freqs, 'Kernel', 'on', 'Location', ...
    'SouthEast', ...
    'Direction', 'out', 'Color', 'b', 'LineStyle', '-', ...
    'LineWidth', 2, 'Marker', '.', 'MarkerSize', 4);

% Set the labels for the axes
xlabel('EUR Allele Frequencies');
ylabel('AFR Allele Frequencies');
title('Allele Frequency Distribution','FontSize',16)

% Edit plot properties.
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'Box', 'off', ...
    'FontWeight', 'normal' ,'XLim',[-0.05,1], 'YLim',[-0.05,1]);

%% Plot the data for all of the genes --- make the plot-

% the colors to use the for two groups
groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10  ] ;

if ~exist('Frequency_Difference_Scatter_10.png','file')
    
    % plot the graph before removing the causal SNPs
    figure()
    hold on
    
    % plot for the EUR
    scatter(snpFreq.af_AFR(snpFreq.highGroup == "EUR"), ...
        snpFreq.af_EUR(snpFreq.highGroup == "EUR"), ...
        5,'filled', ...
        'MarkerFaceColor',groupColors(2,:),'Marker','o',...
        'MarkerEdgeColor',groupColors(2,:),'MarkerFaceAlpha',0.7)
    
    % for the AFR
    scatter(snpFreq.af_AFR(snpFreq.highGroup == "AFR"), ...
        snpFreq.af_EUR(snpFreq.highGroup == "AFR"),...
        5,'filled', ...
        'MarkerFaceColor',groupColors(1,:),'Marker','o',...
        'MarkerEdgeColor',groupColors(1,:))
    
    % plot the non significant data
    nonSigFisher = snpFreq.FDR > 0.05 ;
    scatter(snpFreq.af_AFR(nonSigFisher), ...
        snpFreq.af_EUR(nonSigFisher), ...
        5,'filled','MarkerFaceColor',[ 0.5 0.5 0.5 ],'Marker','o',...
        'MarkerEdgeColor',[ 0.5 0.5 0.5 ],'MarkerFaceAlpha',0.7)
    
    % edit the chart elements
    set(gca,'FontSize',14,'LineWidth',1.5,'Box','off')
    ylabel('SNP Frequency in EUR')
    xlabel('SNP Frequency in AFR');
    legend({sprintf('EUR: %d SNPs',sum(snpFreq_less.highGroup=="EUR")),...
        sprintf('AFR: %d SNPs',sum(snpFreq_less.highGroup=="AFR")),...
        sprintf('Non Sig: %d SNPs',sum(snpFreq_less.highGroup=="None")),...
        },'Location','best')
    title('All Genetic Variants','FontSize',16,'FontWeight','bold')
    
    hold off
    
    % save the figure
    saveas(gca,'Frequency_Difference_Scatter_10.fig','fig')
    saveas(gca,'Frequency_Difference_Scatter_10.png','png')
    
    % print the percentage of variants that differ in frequency
    all_high_eur = sum(snpFreq_less.highGroup=="EUR")/height(snpFreq_less) ;
    all_high_afr = sum(snpFreq_less.highGroup=="AFR")/height(snpFreq_less) ;
    all_none_sig = sum(snpFreq_less.highGroup=="None")/height(snpFreq_less) ;
    
    fprintf('\nPercent higher in EUR variants = %d',all_high_eur*100)
    fprintf('\nPercent higher in AFR variants = %d',all_high_afr*100)
    fprintf('\nPercent none sig variants = %d\n',all_none_sig*100)
    
    fprintf('\nThe number of higher in EUR variants = %d', ...
        sum(snpFreq_less.highGroup=="EUR"))
    fprintf('\nThe number of higher in AFR variants = %d', ...
        sum(snpFreq_less.highGroup=="AFR"))
    fprintf('\nThe number of non sig variants = %d\n', ...
        sum(snpFreq_less.highGroup=="None"))
    
    fprintf('\nAll sig variants btn AFR and EUR = %d\n', ...
        (all_high_eur + all_high_afr)*100)
    
    clear nonSigFisher
end

%% Plot some histogram of variants 

% plot the histogram for europeans
figure()
histogram( log10(snpFreq.af_EUR), 100, 'FaceColor', groupColors(2,:) )
hold on 
% plot the histogram of africans
histogram( log10(snpFreq.af_AFR), 100,'FaceColor', groupColors(1,:))

% edit the plots 
xlabel('Bins (Log 10 of SNP Frequency)')
ylabel('Frequency')
set(gca,'FontSize',14,'LineWidth',2,'Box','off')
legend({'EUR','AFR'},'Location','best')

hold off

% save the figure
saveas(gcf,'Bins_of_SNPs.fig','fig')
saveas(gcf,'Bins_of_SNPs.png','png')

%% Count the instances of Sig Frequencies

% the colors to use the for two groups
groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10 ] ;

% sort the table
snpFreq = sortrows(snpFreq,'freqDiff','descend');

% plot a line graph fo the differences
figure()
% for europeans
bar(snpFreq.af_EUR,'FaceColor',groupColors(2,:))
hold on 
% plot for africans
yyaxis left
bar(-snpFreq.af_AFR,'FaceColor',groupColors(1,:))

% plot the difference
plot(snpFreq.freqDiff,'LineWidth',5,'Color',[0.00,0.45,0.74])

pause(0.5)
set(gca,'YDir','reverse')

% change the figure properties
set(gca,'FontSize',14,'LineWidth',2,'Box','off')
ylabel('Frequency')
xlabel('Number of SNPs')
title('SNP Frequency','FontSize',16,'FontWeight','bold')
legend({'EUR','AFR','Freq Difference'},'Location','best')

hold off

% save the figure
saveas(gcf,'SNP_frequency_line.png','png')

% *********** plot a line graph fo the differences *****************
figure()
% for europeans
bar(snpFreq.af_EUR,'FaceColor',groupColors(2,:))
hold on 
% plot for africans
bar(snpFreq.af_AFR,'FaceColor',groupColors(1,:))
% plot the difference
plot(snpFreq.freqDiff,'LineWidth',5,'Color',[0.00,0.45,0.74])

% change the figure properties
set(gca,'FontSize',14,'LineWidth',2,'Box','off')
ylabel('Frequency')
xlabel('Number of SNPs')
title('SNP Frequency','FontSize',16,'FontWeight','bold')
legend({'EUR','AFR','Freq Difference'},'Location','best')

% save the figure
saveas(gcf,'SNP_frequency_line_2.png','png')

hold off

%% Plot Some Tiled Figures 

% here is the figures
figure(); clf
set(gcf,'position',[500,500,1000,700]);
theLetters = 'a':'c' ;

% make the tile
t = tiledlayout(3,1,'TileSpacing','Compact');

% Tile 1 for the europeans
nexttile
bar(snpFreq.af_EUR,'FaceColor',groupColors(2,:))

% add a letter to the figure and remove the letter from the array
title('SNP Frequency in EUR')
ylabel('Frequency')
% change the figure properties
set(gca,'FontSize',12,'LineWidth',1,'Box','off')
text(-0.07, 1 ,theLetters(1),'Units','normalized','FontWeight','bold', ...
    'FontSize',28)
theLetters(1) = [];

% convert the axis labels
b = get(gca,'XTick')/1e3 ;
b = split( cellstr(num2str(b)) ) ;
b =  strcat(b, 'K')';
set(gca,'XTickLabel',b )

% Tile 2 for the africans
nexttile
bar(snpFreq.af_AFR,'FaceColor',groupColors(1,:))
title('SNP Frequency in AFR')
ylabel('Frequency')

% change the figure properties
set(gca,'FontSize',12,'LineWidth',1,'Box','off')
text(-0.07, 1 ,theLetters(1),'Units','normalized','FontWeight','bold', ...
    'FontSize',28)
theLetters(1) = [];

% convert the axis labels
b = get(gca,'XTick')/1e3 ;
b = split( cellstr(num2str(b)) ) ;
b =  strcat(b, 'K')';
set(gca,'XTickLabel',b )

% Tile 3 for the diffence in snps
nexttile
plot(snpFreq.freqDiff,'LineWidth',5,'Color',[0.00,0.45,0.74])
title('Frequency Difference')
ylabel('Frequency: EUR - AFR')

% change the figure properties
set(gca,'FontSize',12,'LineWidth',1,'Box','off', ...
    'XLim',[-0.2,height(snpFreq)])
text(-0.07, 1 ,theLetters(1),'Units','normalized','FontWeight','bold', ...
    'FontSize',28)
theLetters(1) = [];

% convert the axis labels
b = get(gca,'XTick')/1e3 ;
b = split( cellstr(num2str(b)) ) ;
b =  strcat(b, 'K')';
set(gca,'XTickLabel',b )

% add a common x label
xlabel(t,'Number of SNPs')

% save the figure
saveas(gcf,'SNP_frequency_graph.png','png')

%% Plot the Venn diagram for AFR and EUR SNPs

% plot the venn diagram of SNPs unique in each population
eurSnps = snpFreq( snpFreq.af_AFR == 0 , :) ;
afrSnps = snpFreq( snpFreq.af_EUR == 0 , :) ;
commonSNPs = snpFreq( snpFreq.af_EUR ~= 0 & snpFreq.af_AFR ~= 0,:)  ;

% plot a venn diagram of the snps draw the venn diagram for the SNPs
% that are intersecting
figure()
venn([200,200],75,'EdgeColor','black')
hold on

% add text to the figure: get numbers of intersects to add to the plots
myNumText(3) = height(afrSnps) ;
myNumText(2) = height(commonSNPs);
myNumText(1) = height(eurSnps) ;

% here are text positions and the titles of the plots
textPos = [0.25, 0.47 , 0.70] ;
for jj = 1:length(textPos)
    text(textPos(jj),0.5,num2str(myNumText(jj)),'Units','normalized',...
        'FontWeight','bold','FontSize',16,...
        'HorizontalAlignment','center')
end

% add the titles to the circle
myGroups = {'EUR','AFR'};
textPos = [0.33 , 0.66] ;
for jj = 1:length(textPos)
    text(textPos(jj), 0.98 ,[myGroups{jj}],'Units','normalized', ...
        'FontWeight','bold','FontSize',14,...
        'HorizontalAlignment','center')
end
hold off

% save the figure
saveas(gcf,'Venn_Diagram_of_Genes.fig','fig')

%% Plot the Venn Diagram of the Genes in Each SNP Set

fprintf('\nProcessing the Genes Unique to Each Ethnic Groups\n')

% start multiple parallel instrances
numOfCores = 40;
if multiCoreHighMem == true
    try
    % work on the entire datasets 
    parpool('local',numOfCores)
    catch
    end
end

% check if the file has already been processed 
if ~exist('SNP_genes_location_real.xlsx','file')

    % get the SNPs unique to each group
    EUROnlyGenes = cellstr(snpFreq.HugoSymbol(snpFreq.highGroup == "EUR"));
    AFROnlyGenes = cellstr(snpFreq.HugoSymbol(snpFreq.highGroup == "AFR"));

    % run this in a parrallel for loop because the dataset is large make an
    % array for subseting the data
    subArray = 1:round(length(AFROnlyGenes)/60):length(AFROnlyGenes) ;

    % subset the data
    for ii = 1:length(subArray)
        try
            % for the first instances
            geneStruct.(sprintf('data%d',ii)) =  ...
                unique( AFROnlyGenes(subArray(ii):subArray(ii+1),:)) ;
        catch
            % for the last instances
            geneStruct.(sprintf('data%d',ii)) =  ...
                unique( AFROnlyGenes(subArray(ii):length(AFROnlyGenes),:)) ;
        end
    end

    % preallocate the table
    loopGenes = [] ;

    % run the analysis through the loop
    parfor jj = 1:length(subArray)

        % get the subset of the snpFreq data
        curGenes = geneStruct.(sprintf('data%d',jj)) ;

        % some of the snps are found near many genes so make those into a
        % gene list
        for ii = 1:length(curGenes)

            % print something to the screen
            if rem(ii,100) == 0
                fprintf('\nProcessing gene # %d of file # %d\n',ii, jj)
            end

            % check if the genes have a comma
            if contains(curGenes(ii),',')
                % split the multiple genes and add them to end of the cell
                % array
                loopGenes = [loopGenes; split(curGenes(ii),',')];
            else
                loopGenes = [loopGenes; curGenes(ii)];
            end
        end

    end

    % make the loop genes the actual file
    AFROnlyGenes = unique(loopGenes) ;

    % remove the rows that have the multiple genes  and return only the
    % unique genes
    AFROnlyGenes(contains(AFROnlyGenes,',')) = [];
    AFROnlyGenes(ismember(AFROnlyGenes,'-')) = [] ;

    % ******************************************************************

    % run this in a parrallel for loop because the dataset is large make an
    % array for subseting the data
    subArray = 1:round(length(EUROnlyGenes)/60):length(EUROnlyGenes) ;

    % subset the data
    for ii = 1:length(subArray)
        try
            % for the first instances
            geneStruct.(sprintf('data%d',ii)) =  ...
                EUROnlyGenes(subArray(ii):subArray(ii+1),:) ;
        catch
            % for the last instances
            geneStruct.(sprintf('data%d',ii)) =  ...
                EUROnlyGenes(subArray(ii):length(EUROnlyGenes),:) ;
        end
    end

    % preallocate the table
    loopGenes = [] ;

    % run the analysis through the loop
    parfor jj = 1:length(subArray)

        % get the subset of the snpFreq data
        curGenes = geneStruct.(sprintf('data%d',jj)) ;

        % some of the snps are found near many genes so make those into a
        % gene list
        for ii = 1:length(curGenes)

            % print something to the screen
            if rem(ii,100) == 0
                fprintf('\nProcessing gene # %d of file # %d\n',ii, jj)
            end

            % check if the genes have a comma
            if contains(curGenes(ii),',')
                % split the multiple genes and add them to end of the cell
                % array
                loopGenes = [loopGenes; split(curGenes(ii),',')];
            else
                loopGenes = [loopGenes; curGenes(ii)];
            end
        end

    end

    % make the loop genes the actual file
    EUROnlyGenes = unique(loopGenes) ;

    % remove the rows that have the multiple genes  and return only the
    % unique genes
    EUROnlyGenes(contains(EUROnlyGenes,',')) = [];
    EUROnlyGenes(ismember(EUROnlyGenes,'-')) = [] ;

    % **********************************************************************

    % plot a venn diagram of the snps draw the venn diagram for the SNPs
    % that are intersecting
    figure()
    venn([200,200],75,'EdgeColor','black')
    hold on

    % get the common genes
    commonGenes = intersect(EUROnlyGenes,AFROnlyGenes ) ;

    % get the likely SNPs which are those located in the same gene but they
    % different
    candidateCommonGenes = snpFreq(ismember(snpFreq.HugoSymbol,...
        commonGenes),:) ;

    % add text to the figure: get numbers of intersects to add to the plots
    myNumText(1) = height(AFROnlyGenes) - length(commonGenes) ;
    myNumText(2) = length(commonGenes);
    myNumText(3) = height(EUROnlyGenes) - length(commonGenes) ;

    % here are text positions and the titles of the plots
    textPos = [0.25, 0.47 , 0.70] ;
    for jj = 1:length(textPos)
        text(textPos(jj),0.5,num2str(myNumText(jj)),'Units','normalized',...
            'FontWeight','bold','FontSize',16,...
            'HorizontalAlignment','center')
    end

    % add the titles to the circle
    myGroups = {'AFR','EUR'};
    textPos = [0.33 , 0.66] ;
    for jj = 1:length(textPos)
        text(textPos(jj), 0.98 ,[myGroups{jj}],'Units','normalized', ...
            'FontWeight','bold','FontSize',14,...
            'HorizontalAlignment','center')
    end
    hold off

    % save the figure
    saveas(gcf,'Venn_Diagram_of_Genes.fig','fig')

    % get the exact only genes
    EUROnlyGenes = setdiff(EUROnlyGenes,commonGenes) ;
    AFROnlyGenes = setdiff(AFROnlyGenes,commonGenes) ;

    % convert to to a table
    EUROnlyGenes = array2table(EUROnlyGenes,'VariableNames',...
        {'EUR SNP Genes'}) ;
    AFROnlyGenes = array2table(AFROnlyGenes,'VariableNames',...
        {'AFR SNP Genes'}) ;
    bothGroupsSNPs = array2table(commonGenes,'VariableNames',...
        {'Both Groups SNP Genes'} );

    % also save to the supplementray fiel
    EUROnlyGenes(ismissing(EUROnlyGenes),:) = [] ;

    % save the data to excel for pathway analysis
    writetable(bothGroupsSNPs,'SNP_genes_location_real.xlsx',...
        'Sheet','Location of SNPs','Range','C1')
    writetable(EUROnlyGenes,'SNP_genes_location_real.xlsx',...
        'Sheet','Location of SNPs','Range','A1')
    writetable(AFROnlyGenes,'SNP_genes_location_real.xlsx',...
        'Sheet','Location of SNPs','Range','B1')

    addpath('/required_files/venn')

    % plot the venn diagram of SNPs unique in each population
    eurSnps = snpFreq( snpFreq.af_AFR == 0 , :) ;
    afrSnps = snpFreq( snpFreq.af_EUR == 0 , :) ;
    commonSNPs = snpFreq( snpFreq.af_EUR ~= 0 & snpFreq.af_AFR ~= 0,:)  ;

    % plot a venn diagram of the snps draw the venn diagram for the SNPs
    % that are intersecting
    figure()
    venn([200,200],75,'EdgeColor','black')
    hold on

    % add text to the figure: get numbers of intersects to add to the plots
    myNumText(1) = height(afrSnps) ;
    myNumText(2) = height(commonSNPs);
    myNumText(3) = height(eurSnps) ;

    % here are text positions and the titles of the plots
    textPos = [0.25, 0.47 , 0.70] ;
    for jj = 1:length(textPos)
        text(textPos(jj),0.5,num2str(myNumText(jj)),'Units','normalized',...
            'FontWeight','bold','FontSize',16,...
            'HorizontalAlignment','center')
    end

    % add the titles to the circle
    myGroups = {'AFR','EUR'};
    textPos = [0.33 , 0.66] ;
    for jj = 1:length(textPos)
        text(textPos(jj), 0.98 ,[myGroups{jj}],'Units','normalized', ...
            'FontWeight','bold','FontSize',14,...
            'HorizontalAlignment','center')
    end
    hold off

    % save the figure
    saveas(gcf,'Venn_Diagram_of_SNPS.fig','fig')

    % save the data to excel for pathway analysis
    writetable(eurSnps,'SNP_genes_location_real.xlsx',...
        'Sheet','EUR Only SNPs')
    writetable(afrSnps,'SNP_genes_location_real.xlsx',...
        'Sheet','AFR Only SNPs')
else
    % Load the Euro and AFR only gene and snps
    eurSnps = readtable('SNP_genes_location_real.xlsx',...
        'Sheet','EUR Only SNPs');
    afrSnps = readtable('SNP_genes_location_real.xlsx',...
        'Sheet','AFR Only SNPs') ;
    
    % ************** Return only the High Confidence SNPs ***************
    fprintf('\nReturn only the high frequency variants\n')
    eurSnps =  eurSnps( ismember( eurSnps.Variant, ...
        high_qual_variants.rsid), :);
    afrSnps =  afrSnps( ismember( afrSnps.Variant, ...
        high_qual_variants.rsid), :);
    % *******************************************************************
end

clear jj textPos myNumText commonGenes bothGroupsSNPs geneStruct ...
    ans ii loopGenes subArray t theVars

%% SNPs Comparison for African with UK Biobank

% here is the location of the data procesed by mamana
% /cbio/projects/010/epigenetics_mamana/results/...
% h3africa_annot_v0.0.1_epi_genes.norm_snpeff.tsv

% load the data
h3aFreq = readtable('h3africa_annot_v0.0.1_epi_genes.norm_snpeff.tsv',...
    'FileType','text') ;

% remove the population that should not be in data those should be western
% saharawi and those with very very low samples counts 
otherPopCounts = readtable(...
    'all.v6.african-ancestry-only.meta.extended.pops.tsv', ...
    'FileType','text','Delimiter','tab') ;
otherPopCounts.Properties.VariableNames = ...
    {'Pop','Study','Ancestry','SampleSize'} ;
otherPopCounts = sortrows(otherPopCounts,'SampleSize','ascend') ;

% return instances that have a samples size of more than 30 
otherPopCounts( otherPopCounts.SampleSize < 30, :) = [] ;

% ## remove the population that Nicky suggested we remove ##
otherPopCounts( ismember( otherPopCounts.Study, 'agvp'), :) = [] ;

% remove those populations from the h3aFreq data 
% here are the population and for the two datasets 
h3aPops = lower( replace( h3aFreq.Properties.VariableNames, ...
    {'_MAF','_AF'}, '') );
mamanaPops = [ lower(otherPopCounts.Pop'), ...
    {'x1kg','afr','gbr','fin','eur'} ] ;

% find the location of those populations that have more than 30 instances
% in the h3aAfr data 
locPopsToKeep = ismember(h3aPops, mamanaPops) ;

h3aFreq = [h3aFreq(:,{'GENE','ID','CHROM','POS','REF','ALT','EFFECT'}), ...
    h3aFreq(:, locPopsToKeep) ] ;

% change the first variable name to HugoSymbol
h3aFreq.Properties.VariableNames([1,2]) = ["HugoSymbol","Variant"];
h3aFreq = innerjoin( epigeneticGenes, h3aFreq) ;

%% produce a clustergram of the data 

% find the 1KG_MAF in the data
loc_1kg = find(ismember(h3aFreq.Properties.VariableNames,'x1KG_MAF'),true);

% preallocate the non LD snps
nonLDsnps = [] ;

% check that the file has already been processed
try
    if ~exist('snps_clust_data.txt','file') && multiCoreHighMem == true

        % get the clustering data
        clustData = h3aFreq(:,[3, loc_1kg:end]) ;

        % also get the UKB snp freq data
        ukbFreq = snpFreq(:,{'Variant','af_AFR','af_EUR'}) ;
        ukbFreq.Properties.VariableNames = {'Variant','UKB-AFR','UKB-EUR'} ;

        % merge the two table
        clustData = innerjoin(clustData,ukbFreq) ;

        % ************ remove the snps that are in strong LD *************
        % Load the lds scores of the snps from the GWAS catalogue based
        % on the PIC2 calculations and return only the rows with R-square
        % values > 0.05
        fprintf('\n Loading the PIC2 LD data for the GWAS catalog \n')
        ldData = readtable("GWAScat Small PICS2 2021-09-24.txt");

        fprintf('\n Processing the LD data \n')
        % get only the required variables
        ldData = ldData(:,{'IndexSNP','LinkedSNP','Rsquare'});
        ldData(ldData.Rsquare < 0.6,:) = [] ;
        ldData = ldData(:,{'IndexSNP','LinkedSNP'});

        % get the unique index snps
        curSnps = unique(ldData.IndexSNP) ;

        % make an array for subseting the data
        subArray = 1:round(height(clustData)/60):height(clustData) ;

        % subset the data
        for ii = 1:length(subArray)
            try
                % for the first instances
                snpStruct.(sprintf('data%d',ii)) =  ...
                    snpFreq(subArray(ii):subArray(ii+1), :) ;
            catch
                % for the last instances
                snpStruct.(sprintf('data%d',ii)) =  ...
                    snpFreq(subArray(ii):height(snpFreq), :) ;
            end
        end

        % run the analysis through the loop
        parfor jj = 1:length(subArray)

            % get the subset of the snpFreq data
            curClustData = snpStruct.(sprintf('data%d',jj)) ;

            % remove the snps with an r-square value of < 0.6 LD
            for ii = 1:length(curSnps)

                fprintf('\nProcessing subdata %d of %d for snp %d\n', ...
                    jj, length(subArray), ii)

                % printf some thing to the screen
                if rem(ii,100) == 0
                    fprintf(['\nLD analysis for SNPs #%d of %d in' ...
                        ' loop %d\n'],...
                        jj, ii, length(curSnps))
                end

                % get teh current snps
                theLinkedSnps = ldData.LinkedSNP( ismember( ...
                    ldData.IndexSNP, curSnps{ii})) ;

                % remove the snps in strong LD
                curClustData( ismember( ...
                    curCustData.Variant,theLinkedSnps), :) = [] ;

                % add to the growing table
                nonLDsnps = [nonLDsnps ; curClustData] ;
            end
        end

        % make the non LD data the clust data
        clustData = nonLDsnps ;

        clear nonLDsnps

        writetable(clustData,'snps_clust_data.txt')

        % ***************************************************************
    elseif multiCoreHighMem

        % load the data on the cluster
        clustData = readtable('snps_clust_data.txt') ;
    else
        % get the clustering data
        clustData = h3aFreq(:,[3, loc_1kg:end]) ;

        % also get the UKB snp freq data
        ukbFreq = snpFreq(:,{'Variant','af_AFR','af_EUR'}) ;
        ukbFreq.Properties.VariableNames = {'Variant','UKB-AFR','UKB-EUR'} ;

        % merge the two table
        clustData = innerjoin(clustData,ukbFreq) ;
    end

catch

    % get the clustering data
    clustData = h3aFreq(:,[3, loc_1kg:end]) ;

    % also get the UKB snp freq data
    ukbFreq = snpFreq(:,{'Variant','af_AFR','af_EUR'}) ;
    ukbFreq.Properties.VariableNames = {'Variant','UKB-AFR','UKB-EUR'} ;

    % merge the two table
    clustData = innerjoin(clustData,ukbFreq) ;

end

% change the variable name of the clust data to make them easy for
% processing usign the filterOutGenes functions
clustData.Properties.VariableNames(1) = "HugoSymbol" ;

% remove the rows with more than 20% of the snps = 0
clustData( sum(clustData{:,2:end} == 0,2) > width(clustData)*0.9, :) = [];

% filter out some of the snps
numOfGenes = 5000 ; % 50000
clustData = filterOutGenes(clustData, numOfGenes);

% here are teh colors
if multiCoreHighMem == false
    heatColours.Blues = cbrewer('seq','Blues',10); % Blue
    heatColours.Reds = cbrewer('seq','Reds',10); % Oranges
    heatColours.Greens = cbrewer('seq','Greens',10); % Greens
    heatColours.Purples = cbrewer('seq','Purples',10);
else
    heatColours.Purples = parula ;
end

% produce a clutergram
cgoAll = clustergram(clustData{:,2:end}*100, ...
    'Colormap', heatColours.Purples, ...
    'ColumnLabels', regexprep( ...
    clustData.Properties.VariableNames(2:end) , '_\w+', ''), ...
    'RowLabels', clustData.HugoSymbol ,'Linkage','complete', ...
    'ColumnPDist','cosine','Standardize',1,'Dendrogram',1.2);

% move the dendograms
cm = struct('GroupNumber',{1,2},'Annotation',{'Time1','Time2'},...
    'Color',{[1 1 0],[0.6 0.6 1]});
set(cgoAll,'ColumnGroupMarker',cm)

% Gaone Are SNP frequencies different between Africans living in Africa and
% those of African ancestry not born in Africa (African pops samples from
% Africa versus UK Biobank and Afr American).

% *******************************************************************
% check for snps in strong LD
clustData.Properties.VariableNames(1) = "Variant" ;

clustData = innerjoin( h3aFreq(:, ...
    {'HugoSymbol','GeneClass','Variant'}), clustData );

fprintf('\n The SNP clustering data has been processed\n')

%% NM:2 - Which SNPs are unique to Africa or skewed in their frequencies 

% 2. Which SNPs are unique to Africa or skewed in their frequencies between
% populations? Take all populations and group them into African and
% non-African -select cut-off frequency as definition for present/absent
% (remove admixed pops)? Response: Sure, I will come up with a cut-off
% point for this analysis. The first part is done already, so it is just
% the ones unique to Africa from the H3A dataset

% compare the frequence for variants between AFR in and EUR populations

% get afr data
afrData = h3aFreq(:, [3, loc_1kg:end]);

% remove the ethiopian population which appears twice
try
    afrData.Ethiopia_MAF = [] ;
catch
end

% also get the UKB snp freq data
ukbFreq = snpFreq(:, {'Variant', 'af_AFR', 'af_EUR'});
ukbFreq.Properties.VariableNames = {'Variant', 'UKB-AFR', 'UKB-EUR'};

% merge the two table
chiCompData = innerjoin(afrData, ukbFreq);

% remove the population that should not be in data those should be western
% saharawi and those with very very low samples counts 
otherPopCounts = readtable(...
    'all.v6.african-ancestry-only.meta.extended.pops.tsv', ...
    'FileType','text','Delimiter','tab') ;
otherPopCounts.Properties.VariableNames = ...
    {'Pop','Study','Ancestry','SampleSize'} ;
otherPopCounts = sortrows(otherPopCounts,'SampleSize','ascend') ;

% get the studies that Nicky suggested we do away with
nickyStudies = otherPopCounts.Pop( ismember( ...
    otherPopCounts.Study, 'agvp'));
    
% ## remove the population that Nicky suggested we remove ##
otherPopCounts( ismember( otherPopCounts.Study, 'agvp'), :) = [] ;

% return instances that have a samples size of more than 30 
otherPopCounts( otherPopCounts.SampleSize < 30, :) = [] ;

% return only the african population in data
otherPopCounts( ...
    ~ismember(otherPopCounts.Ancestry, 'African-Ancestry'), :) = [];

% here is the number of samples
ukbAFRsamples = 5978;
ukbEURsamples = 383471;

% add the UKB data to the table
ukbCount = { 'ukb-afr', 'UKB', 'African-Ancestry', ukbAFRsamples; 
             'ukb-eur', 'UKB', 'Non-African-Ancestry', ukbEURsamples };

% add these variable to the mamanaPopCounts
otherPopCounts = [otherPopCounts; cell2table(ukbCount, ...
    'VariableNames', otherPopCounts.Properties.VariableNames)];
otherPopCounts.Pop = lower(otherPopCounts.Pop);

% rename columns
chiCompData.Properties.VariableNames(2:end) = ...
    lower(replace(chiCompData.Properties.VariableNames(2:end),...
    {'_AF', '_MAF'}, ''));

% find the location of those populations that have more than 30 instances
% in the h3aAfr data 
locPopsToKeep = intersect( chiCompData.Properties.VariableNames, ...
    otherPopCounts.Pop,'stable') ;
chiCompData = [ chiCompData(:,{'Variant'}), chiCompData(:,locPopsToKeep) ];
otherPopCounts = otherPopCounts( ismember( ...
    otherPopCounts.Pop , locPopsToKeep ), :) ;

% Assertion to ensure population names in otherPopCounts match chiCompData
% after renaming
pops_in_mamana = unique(otherPopCounts.Pop);
pops_in_chiCompData = chiCompData.Properties.VariableNames(2:end);
assert(all(ismember(pops_in_mamana, pops_in_chiCompData)), ...
    'Mismatch between pop names in mamanaPopCounts and chiCompData.');

% Filter instances where frequency of variants is greater than 0.01 in at
% least one AFR population or the UK population
filteredChiCompData = chiCompData(any(chiCompData{:, 2:end} > 0.01, 2), :);

% Extract the populations
populations = filteredChiCompData.Properties.VariableNames(2:end);

if ~exist('snpFreq_diff_fisher_results.csv','file')
    
    % Initialize the results matrix
    numPops = length(populations);
    resultsMatrix = zeros(numPops);
    
    % Loop through each pair of populations and apply Fisher's Exact Test
    for ii = 1:numPops
        for jj = 1:numPops
            if ii ~= jj
                
                % Initialize counters
                counter = 0;
                
                fprintf(['\nRunning analysis for %s population vs ', ...
                    '#%d of %d vs %s population #%d\n'], ...
                    upper(populations{ii}), ii, numPops, ...
                    upper(populations{jj}) ,jj)
                
                % Loop through each variant
                parfor variant_index = 1:size(filteredChiCompData, 1)
                    
                    % Get the frequencies of the variant in the two
                    % populations
                    freqPop1 = filteredChiCompData{variant_index, ii+1};
                    freqPop2 = filteredChiCompData{variant_index, jj+1};
                    
                    % Calculate number of variants in each population
                    numVariantsPop1 = round(freqPop1 * ...
                        otherPopCounts.SampleSize(strcmp( ...
                        otherPopCounts.Pop, populations{ii})));
                    numVariantsPop2 = round(freqPop2 * ...
                        otherPopCounts.SampleSize(strcmp( ...
                        otherPopCounts.Pop, populations{jj})));
                    
                    % Set up the contingency table
                    tempTable = [numVariantsPop1, ...
                        otherPopCounts.SampleSize(...
                        strcmp(otherPopCounts.Pop, populations{ii})) - ...
                        numVariantsPop1 ;
                        numVariantsPop2, otherPopCounts.SampleSize( ...
                        strcmp(otherPopCounts.Pop, populations{jj})) - ...
                        numVariantsPop2 ] ;
                    
                    % Perform Fisher's exact test
                    [~, pValue] = fishertest(tempTable, 'Tail', 'both');
                    
                    % Check if p-value is significant at alpha = 0.05
                    if pValue < 0.05
                        counter = counter + 1;
                    end
                end
                
                % Save the counter in the results matrix
                resultsMatrix(ii, jj) = counter;
            end
        end
    end
    
    % Create a table from the results matrix
    fisherSnpTable = array2table(resultsMatrix, 'VariableNames', ...
        populations,'RowNames', populations);
    
    % Display the table
    disp(fisherSnpTable);
    
    % save the table
    writetable( fisherSnpTable, 'snpFreq_diff_fisher_results.csv')
    
else
    % load the results 
    fisherSnpTable = readtable('snpFreq_diff_fisher_results.csv') ;
end

% ******** remove the studies that nicky suggested we do away with ******
locToKeep = ~ismember( fisherSnpTable.Properties.VariableNames, ...
    lower(nickyStudies) );
fisherSnpTable = fisherSnpTable(locToKeep, locToKeep ) ;
% ***********************************************************************

% here are the population 
populations = replace( upper(fisherSnpTable.Properties.VariableNames), ...
    '_','-') ;

% plot the heatmap 
figure()
h0 = heatmap( populations, populations, fisherSnpTable{:,:}, ...
    'ColorbarVisible','on','ColorScaling','scaled') ;
h0.title('Variant frequency difference')

% Find the colorbar object
h0.CellLabelFormat = '%.0f';

figure()
% another plot the heatmap with percetanges 
h1 = heatmap( populations, populations, ...
    fisherSnpTable{:,:}./height(filteredChiCompData)*100, ...
    'ColorbarVisible','on','ColorScaling','scaled') ;
h1.title('Percentage of variant frequency differences')

% Find the colorbar object
h1.CellLabelFormat = '%.1f';

% get the median percentage of variants that differ
afr_pops = otherPopCounts.Pop( ...
    ismember( otherPopCounts.Ancestry , 'African-Ancestry') ) ;
afr_pops( ismember( afr_pops, 'UKB-AFR') ) = [] ;

% afr mean and min, and max for the snps 
locAfricaData = ismember(fisherSnpTable.Properties.VariableNames,afr_pops);
fisherSnpAfr = fisherSnpTable{ locAfricaData, locAfricaData};
fisherSnpAfr = fisherSnpAfr./height(filteredChiCompData)*100 ;

% remove the zeroes because those are in the diagnals
fisherSnpAfr = fisherSnpAfr(:) ;
fisherSnpAfr(fisherSnpAfr == 0) = [] ;

% get the values 
afr_median = round( median(fisherSnpAfr) , 3)
afr_max = round( max(fisherSnpAfr), 3 )
afr_min = round( min(fisherSnpAfr), 3 )

% get the mean for the EUR vs AFR group
fisherSnpEUR = fisherSnpTable.ukb_eur./height(filteredChiCompData)*100 ;

% remove the "ukb_afr" population from the comparsion
locUKB_AFR = find( ismember( fisherSnpTable.Properties.VariableNames, ...
    'ukb_afr'),true ) ;

% remove that and remove the zero
fisherSnpEUR(locUKB_AFR) = [];
fisherSnpEUR( fisherSnpEUR  == 0) = [] ;

% get the median
eur_median = round( median(fisherSnpEUR) ,3)
eur_max = round( max(fisherSnpEUR),3 )
eur_min = round( min(fisherSnpEUR),3 )


%% Perform the second but better tSNE analysis 

% get afr data
afrData = h3aFreq(:, [3, loc_1kg:end]);

% remove the ethiopian population which appears twice
try
    afrData.Ethiopia_MAF = [] ;
catch
end

% also get the UKB snp freq data
ukbFreq = snpFreq(:, {'Variant', 'af_AFR', 'af_EUR'});
ukbFreq.Properties.VariableNames = {'Variant', 'UKB-AFR', 'UKB-EUR'};

% merge the two table
tsneData = innerjoin(afrData, ukbFreq);

% remove the population that should not be in data those should be western
% saharawi and those with very very low samples counts 
otherPopCounts = readtable(...
    'all.v6.african-ancestry-only.meta.extended.pops.tsv', ...
    'FileType','text','Delimiter','tab') ;
otherPopCounts.Properties.VariableNames = ...
    {'Pop','Study','Ancestry','SampleSize'} ;
otherPopCounts = sortrows(otherPopCounts,'SampleSize','ascend') ;

% ## remove the population that Nicky suggested we remove ##
otherPopCounts( ismember( otherPopCounts.Study, 'agvp'), :) = [] ;

% return instances that have a samples size of more than 30 
otherPopCounts( otherPopCounts.SampleSize < 30, :) = [] ;

% here is the number of samples
ukbAFRsamples = 5978;
ukbEURsamples = 383471;

% add the UKB data to the table
ukbCount = { 'ukb-afr', 'UKB', 'African-Ancestry', ukbAFRsamples; 
             'ukb-eur', 'UKB', 'Non-African-Ancestry', ukbEURsamples };

% add these variable to the mamanaPopCounts
otherPopCounts = [otherPopCounts; cell2table(ukbCount, ...
    'VariableNames', otherPopCounts.Properties.VariableNames) ];

% also add the 1000G to table and AFR and EUR form the 1000G 
% add the UKB data to the table
x1KGCount = {'1KG-AFR', 'x1KG', 'African-Ancestry', 660; 
             '1KG-EUR','x1KG', 'Non-African-Ancestry', 500;
             'x1KG','x1KG', 'Non-African-Ancestry', 2500};
         
% add these variable to the mamanaPopCounts
otherPopCounts = [otherPopCounts; x1KGCount ] ;        
                 
% convert to lower case           
otherPopCounts.Pop = lower(otherPopCounts.Pop);

% save the results to the supplementary data 
writetable( otherPopCounts, 'Supplementary Data 1.xlsx','Sheet',...
    'Genotype Data Populations') 

% get the tsne data
tsneData = rows2vars(tsneData,'VariableNamesSource','Variant') ;
tsneData.Properties.VariableNames(1) = "Population" ;

% rename the 1000G populations 
tsneData.Population = replace( tsneData.Population , ...
    {'AFR_AF','EUR_AF'}, {'1KG-AFR_AF','1KG-EUR_AF'} ) ;

% conver the groups to categorical
tsneData.Population = categorical(tsneData.Population) ;

% fill the missing variables
tsneData{:,2:end} = fillmissing(tsneData{:,2:end},'nearest') ;

% make the colours of the plots
rng(6)
tsneColours = rand(height(tsneData),3) ;

% here are the distance metrics to try out 
disMetrics = {'euclidean','chebychev',...
    'minkowski','cosine','correlation','spearman','hamming','jaccard',...
    'cityblock','seuclidean'} ;

% creating a tiled layout
figure()
tiledlayout('flow', 'TileSpacing', 'compact');

for ii = 1:length(disMetrics)
    
    rng(1) % for reproducibility
    % let try tnse with various distance metris
    yValues = tsne( tsneData{:,2:end}, 'Distance', disMetrics{ii} ) ; 
    
    % here is the plot
    nexttile
   
    gscatter(yValues(:,1),yValues(:,2), tsneData.Population ,...
        tsneColours,'..',20)
    set(gca,'LineWidth',1,'Box','off','FontWeight','normal')
    xlabel('tSNE-1') ;
    ylabel('tSNE-2') ;
    legend('Off')
    title(disMetrics{ii},'FontSize',14','FontWeight','bold')
end

hold off

% add the populations and save the results to table
tsneResults = table( yValues(:,1), yValues(:,2), ...
    tsneData.Population, 'VariableNames',{'tSNE-1','tSNE-2','Population'});

tsneResults.Population = upper(tsneResults.Population) ;
otherPopCounts.Pop = upper(otherPopCounts.Pop) ;

% clean up the population in myTable 
tsneResults.Population = upper( replace( cellstr( ...
    tsneResults.Population), {'_MAF','_AF'}, '') );

assert(all(ismember(cellstr(tsneResults.Population), ...
    otherPopCounts.Pop)),...
    'Mismatch between pop names in tsneResults and mamanaPopCounts.')

tsneResults = innerjoin(tsneResults, otherPopCounts , ...
    'LeftKey','Population', ...
    'RightKey','Pop')  ;

tsneResults.Ancestry = replace( tsneResults.Ancestry, ...
    {'African-Ancestry','Non-African-Ancestry'},{'African','Others'} ) ;

% save the results to excel 
writetable(tsneResults,'tnse tableau populations.xlsx','Sheet','allPops')

% replot the figure now using only the City block data
figure()
% here is the scatter 
gscatter(tsneResults{:,1},tsneResults{:,2}, tsneResults.Ancestry ,...
    [0.4700 0.6700 0.1900; 0.64,0.08,0.18],'..',40)

% add the labels
% ****************** VECTORISED IMPLEMENTATION *******************
% add number the top of the bar graph
text( tsneResults{:,1}, tsneResults{:,2}, tsneResults.Population, ...
     'FontSize',10,'vert','bottom','horiz','left'); 
% ****************** VECTORISED IMPLEMENTATION *******************

% edit some figure properties 
set(gca,'LineWidth',1.5,'Box','off','FontWeight','normal','FontSize',12)
xlabel('tSNE-1') ;
ylabel('tSNE-2') ;
title('Clustering of Samples','FontSize',14','FontWeight','bold')

% save a copy of the tsne data to use from sig SNP comparison plot 
writetable(tsneData, 'clean_snp_tsne_data.txt') 

clear x1KGCount ukbCount yValues disMetrics tsneColours afrData ...
    tsneData

%% *********************** More Plots ************************

% Check a set of the interesting SNPs discussed -up to about 20 or 30 and
% look at the regional freq differences within different African
% populations. A figure like this would be nice

% load the requried for making the plots 
snp_freq_tsne_data = readtable('clean_snp_tsne_data.txt') ;
ukb_traits = readtable('ukb_traits_epigenetic.txt') ;

% transpose the clean_snp_freq_t_tsne data adn clean the data 
snp_freq_tsne_data = rows2vars( snp_freq_tsne_data , ...
    'VariableNamesSource','Population') ;
snp_freq_tsne_data.Properties.VariableNames(1) ="Variant";

% return only data of africans
snp_freq_tsne_data.Properties.VariableNames = replace( ...
   snp_freq_tsne_data.Properties.VariableNames, {'_MAF','_AF'} ,'') ;
afr_pops = tsneResults.Population( ismember( tsneResults.Ancestry , ...
    'African')) ;
snp_freq_tsne_data = [snp_freq_tsne_data(:,1) , snp_freq_tsne_data( ...
    :, ismember( snp_freq_tsne_data.Properties.VariableNames, afr_pops))];

% remove the rows where all the samples have frequencies less than 0.05
% thus only returning the common variants 
snp_freq_tsne_data( all( snp_freq_tsne_data{:,2:end} < 0.05, 2) ,:) = [];

% get the variance of hte snps and add that to the table 
columnwiseVariance = var(snp_freq_tsne_data{:,2:end},0, 2) ;
snp_freq_tsne_data = addvars( snp_freq_tsne_data, ...
    columnwiseVariance, 'After', 'Variant', ...
    'NewVariableNames','freqVariance') ;

% clean the ukb_traits data 
ukb_traits = removevars(ukb_traits,{'af_AFR','af_EUR','gwasDisease'} ) ;
ukb_traits_snp_freq = innerjoin( ukb_traits, snp_freq_tsne_data) ;

% % get only the significant variants in AFrican 
% ukb_traits_snp_freq = ukb_traits_snp_freq( ...
%     ukb_traits_snp_freq.rep_pvalue < 0.05, :) ;

% return only QTLs in the data by removing the non eqtls in the data
ukb_traits_snp_freq = ukb_traits_snp_freq( any( ismember( ...
    ukb_traits_snp_freq{:, {'eqtl','hqtl','sqtl','mqtl'} },'true' ) ...
    , 2), :) ;

% % sort the variants based on the variance 
% ukb_traits_snp_freq = sortrows( ukb_traits_snp_freq, ...
%    'UKB_description','ascend') ;

%% Now plot the figure

% Extract the frequency data excluding some of the variable names 
loc_freq_variance = find( ismember( ...
    ukb_traits_snp_freq.Properties.VariableNames, 'freqVariance') == 1) ;
pop_start_pos = loc_freq_variance + 1; 

% here the number of rows to plots 
num_of_rows = 40 ;

% **************** most occuring variant in the table *****************

% Use tabulate to get the frequency of each variant
% Define the minimum distance between variants
minDistance = 1000000; % for example, 1,000,000 base pairs

% Use tabulate to get the frequency of each variant
variant_counts = tabulate(ukb_traits_snp_freq.Variant);

% Create a new table with variant names and their counts
most_variants = table(variant_counts(:,1), ...
    cell2mat(variant_counts(:,2)), ...
    'VariableNames', {'Variant', 'Count'});

% Assuming the original data is in a table called ukb_traits_snp_freq and
% most_variants is already created and contains Variant names

% Step 1: Add Chromosome and Position to the most_variants table This is
% done by matching the Variant in most_variants with the Variant in
% ukb_traits_snp_freq
[~, loc] = ismember(most_variants.Variant, ukb_traits_snp_freq.Variant);
most_variants.Chrom = ukb_traits_snp_freq.Chrom(loc);
most_variants.Position = ukb_traits_snp_freq.Position(loc);
most_variants.Position = int64(most_variants.Position);

% Sort the table by count in descending order to have the most occurring
% variants on top
most_variants = sortrows(most_variants, 'Count', 'descend') ;

% Step 3: Filter out variants that are within 1MB of each other on the same
% chromosome
for i = 1:height(most_variants)
    
   % Check if the current index exceeds the table size (after deletions)
    if i > height(most_variants)
        fprintf('\nDone processing all the variants\n');
        break; % Exit the loop
    end
    
    current_variant = most_variants(i, :);
    
    % Find indices of variants on the same chromosome and within the
    % minDistance
    same_chrom = most_variants.Chrom == current_variant.Chrom;
    close_position = abs(most_variants.Position - ...
        current_variant.Position) <= minDistance;
    
    % Set the current position to false
    close_position(i) = false;
    
    % Combine conditions, excluding the current variant
    close_variants = same_chrom & close_position;
    
    % Remove the variants that are too close to the current one
    most_variants(close_variants,:) = [];
end

% Display the filtered table
most_variants = sortrows(most_variants,'Count','descend') ;
% *********************************************************************

% ***************** Select Variants to Plot methods **********************
rng(6) % get some random rows to plot % Random Search Method
randi_rows = randi( height(ukb_traits_snp_freq), num_of_rows,1) ;
freqData = ukb_traits_snp_freq{randi_rows, pop_start_pos:end};

% Assuming most_variants is a table with VariantName as the first column
% and K is the number of top rows you want to select
top_variants = most_variants.Variant(1:num_of_rows);

% Initialize a vector to store the indices of the first occurrence
unique_row_indices = zeros(num_of_rows, 1);

% Loop over each variant to find the index of its first occurrence
for i = 1:num_of_rows
    variant = top_variants{i};
    index = find(strcmp(ukb_traits_snp_freq.Variant, variant), 1, 'first');
    unique_row_indices(i) = index;
end

% Most associated variant methods 
freqData = ukb_traits_snp_freq{ unique_row_indices, pop_start_pos:end};
thePlotSource = ukb_traits_snp_freq(unique_row_indices, :) ;

% ***********************************************************************

% Create a list of SNP identifiers from the 'rep_Variant' column for
% labeling the x-axis
snpLabels = thePlotSource.Variant(1:num_of_rows);

% Create a list of population group names from the table header
populationGroups = thePlotSource.Properties.VariableNames( ...
    pop_start_pos:end);

% using colour brewer
colors = cbrewer('qual','Set1',numel(populationGroups)) ;
colors(colors>1) = 1;

% rng(7) % for reproducibility 
% % Assign colors to each population group
% colors = rand(numel(populationGroups),3); 

% Create the scatter plot
figure()
axes('position',[0.20, 0.15, 0.74, 0.56]);
hold on;

% Initialize handles array for the legend
handles = gobjects(numel(populationGroups), 1);

% Loop through each population group and plot the data
for ii = 1:numel(snpLabels)
    
    % get the yData
    yData = freqData(ii, :);
    xData = repelem(ii, numel(populationGroups));
    
    % Plot each population group with a different color
    for jj = 1:numel(populationGroups)
        
        % Only add to legend if it's the first SNP
        if ii == 1
            handles(jj) = scatter(xData(jj), yData(jj), 80, ...
                colors(jj, :), 'filled');
        else
            scatter(xData(jj), yData(jj), 80, colors(jj, :), 'filled');
        end
        
        % Draw vertical lines at each X point
        line([ii, ii], [min(yData), max(yData)],'Color', [0.3 0.3 0.3],...
            'LineStyle', '-','LineWidth',0.5);
    end
end

% Customize the plot
set(gca, 'XTick', 1:numel(snpLabels), 'XTickLabel', ...
    snpLabels, 'LineWidth', 1, 'FontSize', 12, ...
    'XTickLabelRotation', 60, 'XLim', [0.5, length(snpLabels)+ 0.5]);

xlabel('SNPs');
ylabel('Minor Allele Freq');

% Use handles for legend
[~,objh] = legend(handles, populationGroups, 'Location', 'bestoutside');

%// set font size as desired
objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 80); %// set marker size as desired

% ***************  here are the plotted rows **************************
% thePlotSource variable has the original data 

% get the eqtls in the data and convert them to logical array 
theQTLs = {'eqtl','hqtl','sqtl','mqtl'}  ;
for ii = 1:length(theQTLs)
    % change to logical 
    thePlotSource.(theQTLs{ii}) = strcmp( ...
        thePlotSource.(theQTLs{ii}), 'true')  ; 
end

% Create a new column 'classQTL' with empty strings
thePlotSource = addvars( thePlotSource,  ...
    repmat("", height(thePlotSource), 1), 'after','hqtl', ...
    'NewVariableNames', 'classQTL') ;

% Iterate over each row
for i = 1:height(thePlotSource)
    
    % check the condition 
    trueCount = sum([thePlotSource.eqtl(i), ...
        thePlotSource.mqtl(i), thePlotSource.hqtl(i), ...
        thePlotSource.sqtl(i)]);
    
    % add those to the table 
    if trueCount == 1
        if thePlotSource.eqtl(i)
            thePlotSource.classQTL(i) = "eQTL";
        elseif thePlotSource.mqtl(i)
            thePlotSource.classQTL(i) = "mQTL";
        elseif thePlotSource.hqtl(i)
            thePlotSource.classQTL(i) = "hQTL";
        elseif thePlotSource.sqtl(i)
            thePlotSource.classQTL(i) = "sQTL";
        end
    elseif trueCount > 1
        thePlotSource.classQTL(i) = "multiQTL";
    end
end

% convert to categorical 
thePlotSource.classQTL = categorical( thePlotSource.classQTL, ...
    {'eQTL','mQTL','hQTL','sQTL','multiQTL'}) ;

% plot the heatmap for continous data and barplot for categorical data
classesColors.qtls = [0 0.447 0.741; 0.494 0.184 0.556;0.635 0.078 0.184; ...
    0.3010 0.7450 0.9330 ; 0.95, 0.1 0.1] ;
barData2 = double(thePlotSource.classQTL);

% here are the colours made the eqtls 
curColours = classesColors.qtls(unique(barData2), :) ;
theCategories = categories(thePlotSource.classQTL) ;
curCategories = theCategories(unique(barData2)) ;

% plot the altra bars 
% add a heatmap on top to show the QTL class
axes('position',[0.20, 0.725, 0.67,0.025]);
ultraBars( barData2', curColours )

% add the name of the genes to the left of heatmap
dim = [0.1, 0.715,0.1,0.025];
annotation('textbox',dim,'String','QTL',...
    'FitBoxToText','on','FontSize',10,'EdgeColor','none',...
    'HorizontalAlignment','right','FontWeight','bold',...
    'VerticalAlignment','middle');
hold off

% create a legend for the colour in the main plot
createLegendInternal( 0.80,0.89, curCategories, curColours, ...
    'QTL',[10, 10] , [0.011 ,0.08])

% get the limit of the scatter plot to use those as the limit of the
% bar graph
theXlim = get(gca,'XLim') ;

% plot the bar graph for the number of associated traits 
axes('position',[0.2, 0.760, 0.67, 0.15]);
traitCountBar = most_variants.Count(1:num_of_rows) ;
bar( most_variants.Count(1:num_of_rows) ,'FaceColor', [0.466 0.674 0.188])

% % add the text to the top of the bar graph
% % ****************** VECTORISED IMPLEMENTATION *******************
% % add number the top of the bar graph
% text( (1:length(traitCountBar))+0.37, traitCountBar+1, ...
%     num2str(traitCountBar), ...
%     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
%     'Rotation', 0, 'FontSize',10) ;
% % ****************** VECTORISED IMPLEMENTATION *******************

% edit the figure axis
ylabel('# of Traits')

% set some figure properties and add title ot the figure
set(gca,'Box','off','YColor',[0 0 0], 'XColor',[0 0 0] ,'XLim',theXlim, ...
    'XTickLabel', [])

hold off;

clear columnwiseVariance loc_freq_variance ii jj snpLabels yData ...
    xData colors populationGroups pop_start_pos i 

%% 4. Look at novel African SNPs. Response: I will do this. 

% get SNPs that have more than 0.01 frequency in AFR 
novelSnps = snpFreq( snpFreq.af_AFR > 0.1 & snpFreq.af_EUR < 0.001, :);

% 6. Functional analysis of impact of African SNPs (novel or significantly
% higher MAF) -Found SNPs varying significantly across Afr and Eur and
% within Afr populations -look more deeply at those with potential
% functional effects. 

% Specify the variant IDs
if ~exist('snps_func_impact.csv','file')
    
    % get the current variants
    variant_ids = novelSnps.Variant ;
    
    % Initialize the table to store information for all SNPs
    snps_func_impact = table();
    
    % Loop through all variant IDs
    for ii = 1:length(variant_ids)
        
        % Current variant ID
        variant_id = variant_ids{ii};
        
        fprintf('\nGetting information for %s SNP #%d of %d\n', ...
            variant_id, ii, length(variant_ids) )
        
        % Ensembl VEP API endpoint
        api_endpoint = ['https://rest.ensembl.org/vep/human/id/',...
            variant_id, '?content-type=application/json'];
        
        % Create HTTP options set for GET request
        options = weboptions('ContentType', 'json');
        
        % Perform the API request and save the response
        try
            response = webread(api_endpoint, options) ;
        catch ME
            warning(['Unable to fetch information for variant ',...
                variant_id]);
            continue
        end
        
        % Extract information from response
        snpId = get_field(response, 'id', 'N/A');
        most_severe_consequence = get_field(response, ...
            'most_severe_consequence', 'N/A');
        
        % Extract information from the first transcript_consequence
        if isfield(response, 'transcript_consequences') && ...
                ~isempty(response.transcript_consequences)
            
            % get the transcripts consquences 
            try
                transcript_consequences = ...
                    response.transcript_consequences{1};
            catch
                transcript_consequences = response.transcript_consequences;
            end
            
            % now extract the informations
            biotype = get_field(transcript_consequences, 'biotype', 'N/A');
            amino_acids = get_field(transcript_consequences, ...
                'amino_acids', 'N/A');
            sift_prediction = get_field(transcript_consequences, ...
                'sift_prediction', 'N/A');
            sift_score = get_field(transcript_consequences, ...
                'sift_score', NaN);
            polyphen_prediction = get_field(transcript_consequences, ...
                'polyphen_prediction', 'N/A');
            impact = get_field(transcript_consequences, ...
                'impact', 'N/A');
            polyphen_score = get_field(transcript_consequences,...
                'polyphen_score', NaN);
            
        else
            % for when the values are not present
            biotype = 'N/A';
            amino_acids = 'N/A';
            sift_prediction = 'N/A';
            sift_score = NaN;
            polyphen_prediction = 'N/A';
            impact = 'N/A';
            polyphen_score = NaN;
        end
        
        % Store the information for the current SNP in the table
        snp_info = table({snpId}, {most_severe_consequence}, {biotype},...
            {amino_acids},{sift_prediction}, sift_score, ...
            {polyphen_prediction}, {impact}, polyphen_score, ...
            'VariableNames', {'snpId', 'most_severe_consequence', ...
            'biotype', 'amino_acids', 'sift_prediction', 'sift_score', ...
            'polyphen_prediction', 'impact', 'polyphen_score'});
        
        % Append the information to the snps_info table
        snps_func_impact = [snps_func_impact; snp_info];
    end
    
    % Optionally, you can save the table to a CSV file
    writetable( snps_func_impact, 'snps_func_impact.csv');
    
else
    % load the data
    snps_func_impact = readtable('snps_func_impact.csv');
end

% add the missing annotation from Ensembl to the table
snp_annots = readtable('epigeneticGenes_annotation.tsv','FileType','text');

% get only the required colomns from the data 
snp_annots_less = snp_annots(:, {'rs_dbSNP','Polyphen2_HVAR_score',...
    'SIFT_score','HGVSp_snpEff'}) ;

% Find indices of snps_func_impact that have missing annotations in the
% SIFT and PolyPhen columns
missing_idxs = isnan(snps_func_impact.sift_score) & ...
    isnan(snps_func_impact.polyphen_score);

% Find the snpIds of missing annotations
missing_snpIds = snps_func_impact.snpId(missing_idxs);

% Find indices of these snpIds in snp_annots_less
[~, ia, ib] = intersect(missing_snpIds, snp_annots_less.rs_dbSNP);

% Extract Polyphen2_HVAR_score and split it
polyphen_scores = cellfun(@(x) str2double(strsplit(x, ';')),...
    snp_annots_less.Polyphen2_HVAR_score(ib), 'UniformOutput', false);

% Extract SIFT_score and split it
sift_scores = cellfun(@(x) str2double(strsplit(x, ';')), ...
    snp_annots_less.SIFT_score(ib), 'UniformOutput', false);

% Extract HGVSp_snpEff and split it
hgvsps = cellfun(@(x) strsplit(x, ';'), ...
    snp_annots_less.HGVSp_snpEff(ib), 'UniformOutput', false);

% Store the highest Polyphen2_HVAR_score in snps_func_impact
snps_func_impact.polyphen_score(missing_idxs(ia)) = ...
    cellfun(@max, polyphen_scores);

% Store the highest SIFT_score in snps_func_impact
snps_func_impact.sift_score(missing_idxs(ia)) = ...
    cellfun(@max, sift_scores);

% Store the first instance of HGVSp_snpEff in snps_func_impact
first_hgvsps = cellfun(@(x) x{1}, hgvsps, 'UniformOutput', false);

% Initialize the column if it does not exist
snps_func_impact.HGVSp_snpEff = cell(height(snps_func_impact), 1); 
snps_func_impact.HGVSp_snpEff(missing_idxs(ia)) = first_hgvsps;

snps_func_impact = removevars( snps_func_impact, ...
    {'polyphen_prediction','sift_prediction','amino_acids'}) ;

% join the novel snps to the table 
snps_func_impact = innerjoin(novelSnps(:, {'HugoSymbol','GeneClass', ...
    'Variant','af_AFR','af_EUR'}),  ...
    snps_func_impact ,'LeftKey','Variant','RightKeys','snpId') ;

snps_func_impact = addvars( snps_func_impact, ...
    snps_func_impact.af_AFR./snps_func_impact.af_EUR , ...
    'After','af_EUR','NewVariableNames',{'freqDiv'});
snps_func_impact = sortrows(snps_func_impact, 'freqDiv','descend') ;

% Display the updated snps_func_impact table
head(snps_func_impact)

writetable(snps_func_impact,'novel_snps_func_impact.xlsx')

%% Plot the VEP daata

% Create tiled layout with 3 plots vertically
tiledlayout(2,5,'padding','compact');

% Variables to plot
vars_to_plot = {'most_severe_consequence', 'biotype', 'impact'};
plotNames = {'Predicted Consequences','BioType','Impact'}  ;
theLetters = {'a', 'b', 'c'};

% Loop through the variables to plot
for ii = 1:length(vars_to_plot)
    
    % Get next tile
    if ii == 1
        % here is the next plot
        nexttile([1,5]);
    elseif ii == 2
        nexttile([1,3]);
    else
        nexttile([1,2]);
    end
    
    % Extract data for the current variable and 
    variable_data = snps_func_impact.(vars_to_plot{ii});
    
    % also replace some of the very long unnecessary categories
    variable_data = replace( variable_data, ...
         {'transcribed_processed_pseudogene', ...
         'transcribed_unprocessed_pseudogene',...
         'unprocessed_pseudogene','processed_pseudogene'}, 'pseudogene'); 
    variable_data = replace(variable_data,...
        {'protein_coding_CDS_not_defined'} , 'protein_coding');
    
    % replace the _ with spaces
    variable_data = replace(variable_data,{'_','variant'}, ' ') ;
    
    % Count the occurrences of each category
    [counts, categories] = groupcounts(variable_data);
    
    % sort the categories 
    [counts,sort_order] = sort(counts,'descend') ;
    categoryCounts = categories(sort_order) ;
    
    % Create bar graph
    bar(categorical(categoryCounts,categoryCounts), counts);
    
    % add the number of the top of each bar graph
    text(1:length(categoryCounts), counts, num2str(counts), ...
        'vert','bottom','horiz','center');
    
    % Add title
    title(plotNames{ii});
    
    % add a letter to the figure and remove the letter from the array
    if ii == 1
        text(-0.07, 1 ,theLetters{ii},'Units','normalized',...
            'FontWeight','bold','FontSize',24)
    elseif ii == 2
        text(-0.12, 1 ,theLetters{ii},'Units','normalized',...
            'FontWeight','bold','FontSize',24)
    else
        text(-0.18, 1 ,theLetters{ii},'Units','normalized',...
            'FontWeight','bold','FontSize',24)
    end
 
    % edit the plot 
    set(gca,'FontSize',12,'LineWidth',1,'Box','off')
    
    % Add Y-label
    ylabel('Count');
    
    % Rotate x-labels for better visibility
    xtickangle(45);
end

% get the data save
snps_func_impact_save = removevars( snps_func_impact, {'polyphen_score',...
    'HGVSp_snpEff','sift_score'} ) ;

% write the to the excel file 
writetable(snps_func_impact_save, 'Supplementary Data 1.xlsx',...
    'Sheet','Functional Impact') ;

%% AFR vs EUR Further Analysis 

% Biobank frequencies Look at novel African SNPs Differentiate by African
% subgroups Functional analysis of impact of African SNPs (novel or
% significantly higher MAF) What are the African-specific SNPs with
% potential functional effects on epigenetic-related genes?

% get only the significant snps in africans
afrSigSnps =  snpFreq( snpFreq.highGroup == "AFR" & ...
    snpFreq.FDR < 0.05, : ) ;

% find the 1KG_MAF in the data
loc_1kg = find(ismember(h3aFreq.Properties.VariableNames,'x1KG_MAF'),true);

% join with the other africans snps
afrSigSnps = innerjoin( afrSigSnps, h3aFreq(:,[3, loc_1kg:end]) ) ;

% find the 1KG_MAF in the data
loc_1kg = find(ismember(...
    afrSigSnps.Properties.VariableNames,'x1KG_MAF'),true);

% get only the requried samples
afrSigSnps = afrSigSnps(:, [1,loc_1kg:end]) ;

% remove the europeans from the data will this this code better when I get
% the table from Mamana
afrSigSnps = removevars(afrSigSnps,[2,18,19,20,]);

% change the variable name of the clust data to make them easy for
% processing usign the filterOutGenes functions
afrSigSnps.Properties.VariableNames(1) = "HugoSymbol" ;
afrSigSnps.HugoSymbol = cellstr(afrSigSnps.HugoSymbol);

% remove the rows with more than 20% of the snps = 0
afrSigSnps(sum(afrSigSnps{:,2:end} == 0,2) > ...
    width(afrSigSnps)*0.9,:) = [];

% filter out some of the snps
numOfGenes = 10000 ; % 50000
afrSigSnps = filterOutGenes(afrSigSnps, numOfGenes);

% here are teh colors
heatColours.Oranges = cbrewer('seq','Oranges',10); % Oranges

% produce a clutergram
cgoAfr = clustergram(afrSigSnps{:,2:end}*100, ...
    'Colormap',heatColours.Oranges, ...
    'ColumnLabels',...
    replace(afrSigSnps.Properties.VariableNames(2:end) ,'_','-'), ...
    'RowLabels', afrSigSnps.HugoSymbol ,'Linkage','complete', ...
    'ColumnPDist','cosine','Standardize',1);

%% What are SNPs associated with disease from the GWAS catalog

% traits are associated with SNPs in the epigenetic-related genes?
if ~exist('snpFreq','var')
    % load the processed datasets of only snpfreqs of epigenetic genes
    snpFreq = readtable('imputed_snp_freq.csv') ;

    % add the sig groups for each genes in the snp freq data
    snpFreq.SigCount = snpFreq.FDR < 0.05 ;
end

% get only the significant snps in africans
afrSigSnps = snpFreq( snpFreq.highGroup == "AFR" & snpFreq.FDR < 0.05, :);
eurSigSnps = snpFreq( snpFreq.highGroup == "EUR" & snpFreq.FDR < 0.05, :);

% print something to the screen
fprintf('\n Loading the GWAS catalog data \n')

% consider the SNP in Ld with the snps in integrated datasets
if ~exist('gwas_ld_considered.tsv','file')

    % load the data
    if multiCoreHighMem == false
        fprintf('\n Loading the GWAS catalog data \n')
        gwas = readtable('/required_files/gwasCatalogue.tsv',...
            'FileType','text');
    else
        % loading GWAS catalog data
        fprintf('\n Loading the GWAS catalog data \n')
        gwas = readtable('gwasCatalogue.tsv','FileType','text');
    end

    % change some variaable
    gwas = movevars(gwas,{'SNPS','MAPPED_GENE'},'Before',1) ;
    gwas.Properties.VariableNames(2) = {'HugoSymbol'};

    % get only the gwas data for epigenetic genes
    gwas = gwas( ismember( gwas.HugoSymbol, snpFreq.HugoSymbol), :) ;

    % change some variable names in the gwas catalog data match those in
    % required by the cross referencing function
    toChange = ["CHR_ID","CHR_POS","SNPS"] ;
    newNames = ["Chrom","Position","Variant"] ;
    for ii = 1:length(toChange)
        locSNPs = find(ismember(gwas.Properties.VariableNames, ...
            toChange{ii}),true) ;
        gwas.Properties.VariableNames(locSNPs) = newNames(ii) ;
    end

    % get only the unique rows in the gwas catalog
    gwas = unique(gwas) ;
    gwas( isnan(gwas.Position) | isnan(gwas.Chrom), :) = [];

    % code for processing the GWAS positon from hg38 to hg19
    smallGwas = gwas(:,["Chrom","Position","Variant"]) ;
    smallGwas.Pos2 = smallGwas.Position + 1 ;
    smallGwas( isnan(smallGwas.Position) | isnan(smallGwas.Chrom), :) = [];
    smallGwas = movevars(smallGwas,{'Pos2'},"Before","Variant") ;
    smallGwas.Chrom = strcat({'chr'} , num2str(smallGwas.Chrom) ) ;
    smallGwas.Chrom = regexprep( smallGwas.Chrom,' ','') ;
    writetable(smallGwas(:,1:3),'gwashg19_liftover_epigenetic.txt', ...
        'WriteVariableNames',false,'delimiter', ' ')

    % load the gwas catalog hg19 coodinates
    gwashg19 = readtable('hglft_genome_288f5_a077f0.bed', ...
        'FileType','Text','ReadVariableNames',false) ;
    gwashg19.(3) = [] ;
    gwashg19.Properties.VariableNames = ["Chrom","hg19","hg38","Returned"];
    gwashg19.Position = str2double( extractBetween(gwashg19.hg38,':','-'));
    gwashg19.Returned = [] ;

    % subtract a one because its has been added
    gwashg19.Position = gwashg19.Position - 1;

    % now include the snps in the data
    gwashg19 = innerjoin(smallGwas,gwashg19 ,"keys",{'Chrom','Position'}) ;
    gwashg19 = unique(gwashg19) ;
    size(gwashg19)

    gwashg19.Pos2 = [] ;
    gwashg19.Chrom = str2double( extractAfter(gwashg19.Chrom,'chr') ) ;

    % merge with the gwas catalog orignal data and remove the Position and
    % replace with the hg19 variable with are the new posistios that match
    % those in the other datasets
    gwas = outerjoin(gwas,gwashg19(:,{'Variant','hg19','hg38'}), ...
        "keys","Variant","MergeKeys",true)  ;
    gwas.Position = gwas.hg19 ;
    size(gwas)

    % ******************************************************************
    % *************************** LD TRIMMING **************************

    % Load the LD data
    ldData = readtable('GWAScat Small PICS2 2021-09-24.txt') ;

    fprintf('\n Processing the LD data \n')
    % get only the required variables
    ldData = ldData(:,{'IndexSNP','LinkedSNP','x_hg19_','Rsquare'});
    head(ldData) % for debugging purposes
    ldData(ldData.Rsquare < 0.4, :) = [] ;

    % get only the data that exist for africans and europeans for the
    % complete datasets and the causal variant dataset just in case
    % something is being missing 
    ldData = ldData( ismember(ldData.IndexSNP, snpFreq.Variant) | ...
        ismember(ldData.LinkedSNP, snpFreq.Variant),:) ;

    % trim the gwas catalogue to only return the variants associated with
    % pulmnarory function - change the snps in the GWAS catelogue data to
    % the ones reported in this study if those snps are in strong LD - all
    % the data is coming from the same file so I can change those here

    % add the chrom and position to the intData
    epiVariantqc = readtable('variantqc_epi.txt');
    epiVariantqc = unique(epiVariantqc) ;
    epiVariantqc.HugoSymbol = [] ;
    epiVariantqc.Properties.VariableNames(1:2) = {'Chrom','Position'} ;
    snpFreq_pos = innerjoin(snpFreq, epiVariantqc, "LeftKey","Variant", ...
        "RightKey","rsid") ;


    [gwasC_Ld, gwas] = epigenetics_ld_crossref(ldData, snpFreq_pos, ...
        gwas, 0.6, 500000) ;

    % *************************** LD TRIMMING **************************
    % ******************************************************************

    fprintf('\n The number of LD variants in GWAS are %d\n', ...
        length(height(gwasC_Ld)))

    % save the results to a file 
    writetable(gwas,'gwas_ld_considered.tsv','FileType','text',...
        'delimiter','tab')

    clear gwashg19 snpFreq_pos
    
else
    % load the processed data 
    gwas = readtable('gwas_ld_considered.tsv','FileType','text');
end

% remove the adjusted for values in the data 
gwas.DISEASE_TRAIT = replace(gwas.DISEASE_TRAIT, ...
    {'adjusted for body mass index','adjusted for BMI'},'');

% cavert the disease traits to categorical
gwas.DISEASE_TRAIT = categorical(gwas.DISEASE_TRAIT) ;

% return only the epigenetic genes
gwas(~ismember(gwas.Variant,snpFreq.Variant), :) = [] ;

% get only the disease and variant from the table
gwasDisease = gwas(:, {'Variant','DISEASE_TRAIT'} ) ;
gwasDisease = unique(gwasDisease ) ;

% add the gwas catalog phenotypes to the snps
gwasAFR = innerjoin( afrSigSnps, gwasDisease,'Key','Variant') ;
gwasEUR = innerjoin( eurSigSnps, gwasDisease,'Key','Variant') ;

% convert the disease traits to categorical and get the disease associated
% with the snps
fprintf('\n Getting the diseases associated with SNPs in each group \n')

afrDisease = table( categories(gwasAFR.DISEASE_TRAIT) , ...
    countcats(gwasAFR.DISEASE_TRAIT) , 'VariableNames', ...
    {'Disease','Afr_Count'}) ;
afrDisease = sortrows(afrDisease,'Afr_Count','descend')  ;

% convert the disease traits to categorical
eurDisease = table( categories(gwasEUR.DISEASE_TRAIT) , ...
    countcats(gwasEUR.DISEASE_TRAIT) , 'VariableNames', ...
    {'Disease','Eur_Count'}) ;
eurDisease = sortrows(eurDisease,'Eur_Count','descend')  ;

% let see all traits in the gwas catalog
gwasDisease = table( categories( gwasDisease.DISEASE_TRAIT), ...
    countcats( gwasDisease.DISEASE_TRAIT) , 'VariableNames',  ...
      {'Disease','All_Count'}) ; 
gwasDisease = sortrows(gwasDisease,'All_Count','descend')  ;

% put the three tables together
gwasDisease = outerjoin( outerjoin( gwasDisease,afrDisease , ...
    'MergeKey',true), eurDisease ,'MergeKey',true) ;
gwasDisease = sortrows(gwasDisease,'All_Count','descend')  ;

% get the top20 disease of each group
top10diseases = [afrDisease.Disease(1:9), eurDisease.Disease(1:9)] ;
top10diseases = gwasDisease(ismember(gwasDisease.Disease,top10diseases),:);

% ***** Which among the gene class is more assoicatd with disease ******

% get the most frequent affected genes class in africans
afr_geneRank = table( categories(gwasAFR.GeneClass), ...
    countcats( gwasAFR.GeneClass) , 'VariableNames',  ...
      {'GeneClass','Afr_Count'}) ; 
afr_geneRank = sortrows(afr_geneRank,'Afr_Count','descend')  ;

% do the same for europeans
eur_geneRank = table( categories(gwasEUR.GeneClass), ...
    countcats( gwasEUR.GeneClass) , 'VariableNames',  ...
      {'GeneClass','Eur_Count'}) ; 
eur_geneRank = sortrows(eur_geneRank,'Eur_Count','descend')  ;

% get the frequency of genes in the each class
all_geneRank = table( categories(epigeneticGenes.GeneClass), ...
    countcats( epigeneticGenes.GeneClass) , 'VariableNames',  ...
      {'GeneClass','All_Count'}) ; 
all_geneRank = sortrows(all_geneRank,'All_Count','descend')  ;

% put the three tables together
all_geneRank = outerjoin( outerjoin( all_geneRank, afr_geneRank, ...
    'MergeKey',true), eur_geneRank ,'MergeKey',true) ;
all_geneRank = sortrows(all_geneRank,'All_Count','descend')  ;

% BROMO and HDM snps are associated with more trait in european than would
% be expected by chance

% print something to the screen 
fprintf('\nthe number of epigenetic variants is %d\n', ...
    length(unique(gwas.Variant)) ) 
fprintf('\nthe number of trait associations is %d\n', ...
    length(unique(gwas.DISEASE_TRAIT)) ) 

% save the results to excel for analysis in tableau
writetable(gwas,'GWAS Diseases.xlsx','Sheet','all_epigenetic_traits') 
writetable(gwasAFR,'GWAS Diseases.xlsx','Sheet','GWAS AFR Diseases') 
writetable(gwasEUR,'GWAS Diseases.xlsx','Sheet','GWAS EUR Diseases') 
writetable(gwasDisease,'GWAS Diseases.xlsx','Sheet','All Diseases Count') 
writetable(top10diseases,'GWAS Diseases.xlsx','Sheet','Top 10 Diseases') 
writetable(afrDisease,'GWAS Diseases.xlsx','Sheet','AFR Disease Count') 
writetable(eurDisease,'GWAS Diseases.xlsx','Sheet','EUR Disease Count') 
writetable(all_geneRank,'GWAS Diseases.xlsx','Sheet','Gene Class Count')

% clear all_geneRank top10diseases gwasDisease afrDisease ...
%     eurDisease all_geneRank

%% Plot some GWAS figures 

% plot a figure of disease association counts in the GWAS catalog for
% epigenetic genes
top10diseases.Disease = regexprep(top10diseases.Disease, ...
    'Post bronchodilator ','');
top10diseases.Disease = regexprep(top10diseases.Disease, ...
    'Sex hormone-binding globulin levels','Sex hormone-binding globulin');
x = categorical(top10diseases.Disease,top10diseases.Disease);
vals = top10diseases{:,2:end};

% the group colours are flipied here 
theLetters = 'a':'j';

% get only teh top 10 values 
x = categorical( cellstr( x(1:10) ), cellstr(x(1:10)) );
vals = vals( 1:10,:);

% plot the bar graph 
figure()

% plot the bar graphs using a tiled layout
tiledlayout(2,1,'padding','compact');

% Plot the data on a bargraph of on the tile
nexttile

% here is the bar graph
b = bar( x, vals);

% edit some plot features 
set(gca,'LineWidth',1,'Box','off','FontWeight','bold','FontSize',12)
ylabel('Trait count') ;
title('Top GWAS Traits Associated with Epigenetic Variants',...
    'FontSize',14','FontWeight','bold')

% add a letter to the figure and remove the letter from the array
text(-0.07, 1 ,theLetters(1),'Units','normalized','FontWeight','bold', ...
    'FontSize',28)

% add numbers to the top of the plots 
for ii = 1:3
    xtips1 = b(ii).XEndPoints;
    ytips1 = b(ii).YEndPoints;
    labels1 = string(b(ii).YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize',12)
    
    % also change the face color
    if ii == 1
        b(ii).FaceColor = [0.7 0.7 0.7] ;
    else
        b(ii).FaceColor = groupColors(ii-1,:) ;
    end
end

% add a legend and title to the figure
legend({'All','Higher in AFR','Higher in EUR'},'Box','off')

% ******************* produce a plot for the genes ********************
if iscell(gwas.HugoSymbol)
    gwas.HugoSymbol = categorical(gwas.HugoSymbol) ;
end

% here are the unique gwas studies 
[~,theUnique] = unique( gwas(:,{'Variant','DISEASE_TRAIT'} ) );
uniqueGwas = gwas(theUnique,:);

% here are all the genes 
gwasGenes = table( categories( uniqueGwas.HugoSymbol), ...
    countcats(uniqueGwas.HugoSymbol),'VariableNames', ...
    {'HugoSymbol','All_Count'});
gwasGenes = sortrows(gwasGenes,'All_Count','descend')  ;

% here are the african genes 
topAFRgenes = table( categories(gwasAFR.HugoSymbol), ...
    countcats(gwasAFR.HugoSymbol), ...
    'VariableNames',{'HugoSymbol','Afr_Count'}) ; 
topAFRgenes = sortrows(topAFRgenes,'Afr_Count','descend')  ;

% here are the whites genes 
topEURgenes = table( categories(gwasEUR.HugoSymbol), ...
    countcats( gwasEUR.HugoSymbol),'VariableNames',  ...
      {'HugoSymbol','Eur_Count'}) ; 
topEURgenes = sortrows(topEURgenes,'Eur_Count','descend')  ;

% put the data together 
% put the three tables together
gwasGenes = outerjoin( outerjoin( gwasGenes,topAFRgenes , ...
    'MergeKey',true), topEURgenes ,'MergeKey',true) ;
gwasGenes = sortrows(gwasGenes,'All_Count','descend')  ;

% here are the top genes
top10genes = [topAFRgenes.HugoSymbol(1:9),topEURgenes.HugoSymbol(1:9)] ;
top10genes = gwasGenes(ismember(gwasGenes.HugoSymbol,top10genes), :) ;

% let see if the genes are a cell array 
if iscell(top10genes.HugoSymbol)
    top10genes.HugoSymbol = categorical(top10genes.HugoSymbol, ...
        top10genes.HugoSymbol) ;
end

% Plot the data on a bargraph of on the tile
nexttile

% plot the bar graph
topGenePlot = top10genes(1:10,:);
topGenePlot.HugoSymbol = categorical(cellstr(topGenePlot.HugoSymbol), ...
    cellstr(topGenePlot.HugoSymbol) );
b = bar( topGenePlot.HugoSymbol, topGenePlot{:,2:end} );

% edit some plot features 
set(gca,'LineWidth',1,'Box','off','FontWeight','bold','FontSize',12)
xlabel('Genes') ;
ylabel('Sum of disease variants') ;
title('\bfTop Epigenetic Genes Associated with GWAS Traits','FontSize',14)

% add numbers to the top of the plots 
for ii = 1:3
    xtips1 = b(ii).XEndPoints;
    ytips1 = b(ii).YEndPoints;
    labels1 = string(b(ii).YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize',12)
    
    % also change the face color
    if ii == 1
        b(ii).FaceColor = [0.7 0.7 0.7] ;
    else
        b(ii).FaceColor = groupColors(ii-1,:) ;
    end
end

% add a letter to the figure and remove the letter from the array
text(-0.07, 1 ,theLetters(2),'Units','normalized','FontWeight','bold', ...
    'FontSize',28)

% add a legend and title to the figure
legend({'All','Higher in AFR','Higher in EUR'},'Box','off')

clear ii x vals  xtips1  ytips1 labels1 gwasGenes topEURgenes ...
    topAFRgenes theUnique topGenePlot 

%% Here are the predictions of SNP pathogenicity 

% https://sites.google.com/site/jpopgen/dbNSFP
% dbNSFP is a database developed for functional prediction and annotation
% of all potential non-synonymous single-nucleotide variants (nsSNVs) in
% the human genome. Its current version is based on the Gencode release 29
% / Ensembl version 94 and includes a total of 84,013,490 nsSNVs and ssSNVs
% (splicing-site SNVs).  It compiles prediction scores from 38 prediction
% algorithms (SIFT, SIFT4G, Polyphen2-HDIV, Polyphen2-HVAR, LRT,
% MutationTaster2, MutationAssessor, FATHMM, MetaSVM, MetaLR, MetaRNN,
% CADD, CADD_hg19, VEST4, PROVEAN, FATHMM-MKL coding, FATHMM-XF coding,
% fitCons x 4, LINSIGHT, DANN, GenoCanyon, Eigen, Eigen-PC, M-CAP, REVEL,
% MutPred, MVP, MPC, PrimateAI, GEOGEN2, BayesDel_addAF, BayesDel_noAF,
% ClinPred, LIST-S2, ALoFT), 9 conservation scores (PhyloP x 3, phastCons x
% 3, GERP++, SiPhy and bStatistic) and other related information including
% allele frequencies observed in the 1000 Genomes Project phase 3 data,
% UK10K cohorts data, ExAC consortium data, gnomAD data and the NHLBI Exome
% Sequencing Project ESP6500 data, various gene IDs from different
% databases, functional descriptions of genes, gene expression and gene
% interaction information, etc.

% Here is the data processed using the snps_annoations_processing.m
% function

% load the epigenetic genes information
if ~exist('epigeneticGenes_annotation.tsv','file')
    
    % get the chromosome data 
    chroms = readtable('chromosome_name.csv','Format','auto', ...
        'ReadVariableNames',false) ;
    chroms = chroms.Var1(1:22) ;

    % preallocate the all data table
    snp_annots = [] ;

    % process the data from chromosome 1 to 23
    for ii = 1:length(chroms)

        % get the current chromosome name
        curChrom = sprintf('epigenes_chr%s.txt', chroms{ii} ) ;

        fprintf('\nLoading data for chromosome %d: %s\n',ii,curChrom)

        % load the current data of the chromosome 
        curData = readtable(['only_epigenes/',curChrom],'FileType','text');
        
        % some files seem to be corrupt so I skip the head rows for those
        if width(curData) < 500
            curData = readtable(['only_epigenes/',curChrom],...
                'FileType','text','HeaderLines',13);
        end
        
        % change the first variable name to HugoSymbol to match the
        % epigentic genes table
        curData.Properties.VariableNames(13) = "HugoSymbol";

        % let see a small part of the current data
        curData(1:5,1:12)

        % find the position of the 'Geuvadis_eQTL_target_gene' var
        locEnd = find( contains( curData.Properties.VariableNames, ...
            'Geuvadis_eQTL_target_gene','IgnoreCase',true),true) ;
        curData = curData(:,1:locEnd);

        % end the loop when I get the x chromosome because the UKB data
        % does not have those datasets
        if iscell(curData.x_chr)
            break
        end

        % also return only the snps with Ids in the dbsnp
        curData.rs_dbSNP = replace(curData.rs_dbSNP,'.','') ;
        curData( cellfun(@isempty,curData.rs_dbSNP), :) = [];
        
        % remove the nan values
        curData(isnan(curData.x_chr),:) = [];

        % remove the unrequired variables
        theVars = curData.Properties.VariableNames ;
        curData = curData(:, ...
            ~contains(theVars,{ ...
            'gNomad','ExAC','x1000','Ensembl','hg18_','codonpos',...
            'clinvar_MedGen_id','clinvar_Orphanet_id','UK10',...
            'clinvar_var_source','clinvar_hgvs','clinvar_id',...
            'codon_degeneracy','AltaiNeandertal','Denisova',...
            'ChagyrskayaNeandertal','VindijiaNeandertal','cds_strand',...
            'VEP_canonical','TSL','GENCODE_basic','APPRIS',...
            'Uniprot_entry','Uniprot_acc'}, ...
            'IgnoreCase',true)) ;

        % convert the numbers to cell array
        for jj = 7:width(curData)
            if isnumeric( curData.(jj) )
                curData.(jj) = cellstr( num2str(curData.(jj))) ;
            end
        end

        % make sure the variable names are the same
        if ii > 1
            if width(snp_annots) ~= width(curData)
                curData = curData(:, ...
                    ismember(curData.Properties.VariableNames,...
                    snp_annots.Properties.VariabelNames)) ;
            end
        end

        % add to growing the data
        snp_annots = [snp_annots; curData] ;

    end

    % remove the columns with 50% missing data
    snp_annots(:, sum(ismissing(snp_annots))/height(snp_annots) > 0.2) =[];

    % let see how much data we have
    fprintf('\nThe complete data have %d rows and %d columns\n', ...
        height(snp_annots) , width(snp_annots) )

    % save the data to text file
    writetable(snp_annots,'epigeneticGenes_annotation.tsv')
else
    % load the processed data
    snp_annots = readtable('epigeneticGenes_annotation.tsv', ...
        'FileType','text') ;
    
    % return only the variants present in snpFreq 
    snp_annots = snp_annots( ismember( snp_annots.rs_dbSNP , ...
        snpFreq.Variant), :) ;
    
end

% merge the data with the snpFreq information - snpFreq

%% ************* Plot some pie charts of the snps ****************** 

pieLabels = snp_annots(:, {'SIFT4G_pred','SIFT4G_converted_rankscore'} );

% the pie labels have multiple labels get those which have atleast 1
% predicted pathogenic variant 
pieLabels.SIFT4G_pred(contains(pieLabels.SIFT4G_pred,'D')) = {'D'} ;
pieLabels.SIFT4G_pred(~contains(pieLabels.SIFT4G_pred,'D')) = {'T'} ;

% convert to categorical
pieLabels.SIFT4G_pred = categorical(pieLabels.SIFT4G_pred) ;
pieLabels.SIFT4G_pred = renamecats(pieLabels.SIFT4G_pred,{'D','T'}, ...
    {'Deleterious','Tolerated'}) ;

% get the goup statistics of the data 
pieLabels = grpstats(pieLabels,'SIFT4G_pred') ;

% sort the pieLabels so that the pie chart can be plotted
% properly
pieLabels = sortrows(pieLabels, 'GroupCount','descend') ;

% create the string variable for the pie chart labels
finalLabels = strcat(cellstr(pieLabels.SIFT4G_pred),...
    strcat(':XX',num2str(pieLabels.GroupCount) ) ) ;

finalLabels = strrep( strrep(finalLabels,' ', '') ,'XX',' ') ;

% specificy the parts of the pie chart to explode
explode = pieLabels.GroupCount / sum(pieLabels.GroupCount) ...
    < 0.05 ;

% finally produce the pie chart
figure()
pie(pieLabels.GroupCount ,explode, finalLabels)

clear ii curMuts cancerStudies

%% How do the eQTLs vary between the tissue in GTEx tissues

% There are difference in which genes are eQTLs depending on the tissues.
% Let's check for the this explicitly

% check if the files have been processed
if ~exist('gtex_epigenetic.xlsx','file')

    % define the tissue types available that are aviable in the folder
    if exist('/required_files/gtexData','dir')
        addpath( ...
            '/required_files/gtexData/GTEx_Analysis_v8_eQTL')
        cd('/required_files/gtexData/GTEx_Analysis_v8_eQTL')
        listing = ...
            dir( ...
            '/required_files/gtexData/GTEx_Analysis_v8_eQTL');
        listing = struct2table(listing) ;
        listing = listing.name(contains(listing.name,'egenes')) ;
    end
    
    % read one Gtex tissue one at a time to one table of signficant eQTLs
    
    % get the tissue names
    tissues = extractBefore(listing,'.') ;
    for ii = 1:length(tissues)
        % add the path
        if ii == 1
            addpath("ExtactedGTExFiles/",'-end')
        end
        
        % checck that the file has been unzipped previously
        if ~exist([tissues{ii} ,'.v8.egenes.txt'], 'file')
            gunzip( listing{ii} ,'ExtactedGTExFiles')
        end
    end
    
    % set the values out side the parfor loop
    eQTLsAll = [] ;
    
    % run the parfor loop
    for ii = 1:length(tissues)
        
        % load the GTEx data of eQTLs
        eqtlsCur = readtable([tissues{ii} ,'.v8.egenes.txt']) ;
        
        % return only the genes involved in epigenetics
        eqtlsCur(~ismember(eqtlsCur.gene_name, ...
            epigeneticGenes.HugoSymbol),:)=[];
        
        % replace the snp that are longer than 2 with X
        locThem = cellfun(@length, eqtlsCur.ref) > 1 ;
        longerSNPs = repmat({'X'}, length(eqtlsCur.ref(locThem)),1 ) ;
        eqtlsCur.ref(locThem) = longerSNPs ;
        
        % also do the same the alternate SNPs
        locThem = cellfun(@length, eqtlsCur.alt) > 1 ;
        longerSNPs = repmat({'X'}, length(eqtlsCur.alt(locThem)),1 ) ;
        eqtlsCur.alt(locThem) = longerSNPs ;
        
        % add the column for the type of alteration
        eqtlsCur.alterationType = strcat(eqtlsCur.ref , ...
            strcat('->',eqtlsCur.alt)) ;
        
        % add the tissue names to the eQTLs
        eqtlsCur = addvars( eqtlsCur , ...
            repmat( tissues(ii), height(eqtlsCur) , 1 ), ...
            'Before',1 , 'NewVariableNames', {'Tissue'} ) ;
        
        % get all eqtls and put them in one table
        eQTLsAll = [eQTLsAll ; eqtlsCur ];
        
        % return only the genes and the alterationType and p values
        eqtlsCur = eqtlsCur(:, ...
            {'gene_name','alterationType','pval_beta'}) ;
        
        % save a copy of all the SNP that are present in the population
        SNPcur = eqtlsCur(:, {'gene_name','alterationType'}) ;
        SNPcur.Properties.VariableNames(1:2) = [ ...
            'HugoSymbol', tissues(ii) ] ;
        
        % also delete the SNPs that are not signicant between the
        % comparisons
        eqtlsCur.alterationType(eqtlsCur.pval_beta > 0.05) = {''} ;
        
        % replace the alteration that are longer with others
        eqtlsCur.alterationType(contains(eqtlsCur.alterationType,'X'))...
            = {'others'} ;
        eqtlsCur.alterationType = categorical(eqtlsCur.alterationType);
        
        % replace the most infrequent mutation with multiple output like
        % the TCGA does in the publication
        
        % return only the first two columns of the table and change the
        % name of the column to tissue name
        eqtlsCur = eqtlsCur(:, {'gene_name','alterationType'}) ;
        eqtlsCur.Properties.VariableNames(1:2) = ...
            ['HugoSymbol', tissues(ii)];
        
        if ii == 1
            gtexEqtls = eqtlsCur ;
            gtexSnps = SNPcur;
        else
            % add to the growing tables
            gtexEqtls = outerjoin(gtexEqtls,eqtlsCur,"MergeKeys",true) ;
            gtexSnps  = outerjoin(gtexSnps,SNPcur,"MergeKeys",true) ;
        end
    end
    % save the files
    writetable(gtexEqtls,'gtex_epigenetic.xlsx','Sheet','eqtls')
    writetable(gtexSnps,'gtex_epigenetic.xlsx','Sheet','snps')

    cd /required_files/epigenetics
else
    % load the processed data
    gtexEqtls = readtable('gtex_epigenetic.xlsx','Sheet','eqtls') ;
    gtexSnps = readtable('gtex_epigenetic.xlsx','Sheet','snps') ;
end

clear longerSNPs SNPcur listing ii eqtlsCur locThem

%% Plot the eQTLs in various tissues to see differences

% delete genes eqtls that have all missing data
gtexSnps( all( ismissing(gtexEqtls{:,2:end}), 2),:) = [] ;
gtexEqtls(all( ismissing(gtexEqtls{:,2:end}), 2),:) = [] ;

% get the data for producing a clustergram of tissues that share eqtls
clustData = double( double( categorical( gtexEqtls{:, 2:end} ) ) > 0 );

% produce a clutergram
cgo = clustergram(clustData','Colormap',redbluecmap,'RowLabels',...
    replace(gtexEqtls.Properties.VariableNames(2:end) ,'_','-' ), ...
    'ColumnLabels', gtexEqtls.HugoSymbol ,'Linkage',...
    'complete','ColumnPDist','euclidean','Dendrogram', 6);

% Set properties of the clustergram object
cgo.addTitle('eQTLs Across Various Tissues');
addXLabel(cgo,'eQTL SNPs','FontSize',12);
addYLabel(cgo,'Tissues','FontSize',12);

%% Process the eQTL data with only p-values for SNP that are eqtls

% check that the file exist
if ~exist('snpWithEqtlsPvalues.xlsx','file')
    % run the parfor loop
    for ii = 1:length(tissues)
        
        % load the GTEx data of eQTLs
        eqtlsCur = readtable([tissues{ii} ,'.v8.egenes.txt']) ;
        
        % return only the genes involved in epigenetics
        eqtlsCur(~ismember(eqtlsCur.gene_name, ...
            epigeneticGenes.HugoSymbol),:)=[];
        
        % return only the genes and the rsID and p values
        eqtlsCur = eqtlsCur(:, ...
            {'gene_name','rs_id_dbSNP151_GRCh38p7','pval_beta'}) ;
        
        % return only the signficant snp ids in the tissues
        eqtlsCur(eqtlsCur.pval_beta > 0.05,:) = []  ;
        
        % save a copy of all the SNP that are present in the population
        eqtlsCur.Properties.VariableNames(1:3) = [ ...
            'HugoSymbol','rsID', tissues(ii) ] ;
        
        % replace the most infrequent mutation with multiple output like
        % the TCGA does in the publication
        if ii == 1
            snp_eqtls_pvalues = eqtlsCur ;
        else
            % add to the growing tables
            snp_eqtls_pvalues = outerjoin( ...
                snp_eqtls_pvalues,eqtlsCur,"MergeKeys",true) ;
        end
    end
    
    % Plot the SNPs that have signifcant eQTLs across various tissues
    
    % delete the columns that have all NaN values
    snp_eqtls_pvalues(:, all( isnan(snp_eqtls_pvalues{:,3:end}),1)) = [];
    
    % save the file to excel
    writetable(snp_eqtls_pvalues,'snpWithEqtlsPvalues.xlsx')
    
else % or just load the processed data
    snp_eqtls_pvalues = readtable('snpWithEqtlsPvalues.xlsx');
end

% make the variable names the same as those in the snps table
snp_eqtls_pvalues.Properties.VariableNames(2) = "Variant" ;

% add the frequency of snps to the table
snp_eqtls_pvalues = innerjoin(snpFreq,...
    removevars(snp_eqtls_pvalues,'HugoSymbol'), 'Key','Variant');

% remove the cells from the data
snp_eqtls_pvalues = removevars(snp_eqtls_pvalues, ...
    {'Cells_Cultured_fibroblasts','Cells_EBV_transformed_lymphocytes', ...
    'Minor_Salivary_Gland'}) ; 

% replace some variable names
snp_eqtls_pvalues.Properties.VariableNames = replace( ...
    snp_eqtls_pvalues.Properties.VariableNames, {'Whole_Blood'}, ...
    'WholeBlood' ) ;

% get the data for producing a clustergram of tissues SNPs that are
% significant eQTLs across various tissues and replace the NaN values with
% 1 to mean not significant although it is not the same thing. This should
% only be donw for the eqlts data remove the variables that are not
% required but also add the mean values of each to the table

% get the start of the eqtls agains
eqtl_start = find( ismember(snp_eqtls_pvalues.Properties.VariableNames,...
    'Adipose_Subcutaneous'),true) ;

% here are the variables to use
theVars = snp_eqtls_pvalues.Properties.VariableNames(eqtl_start:end) ;

% loop over the vars
for ii = 1:length(theVars) 
    
    % get the current variable names
    curVar = theVars(ii) ;
    
    % check if the variable has a _
    if contains(curVar,'_')
        curVar = extractBefore(curVar,'_') ;
    else
        continue
    end
    
    % get all the variables that contains that particular snps
    locVars = contains(snp_eqtls_pvalues.Properties.VariableNames,curVar);
    cur_pvalues = min( snp_eqtls_pvalues{:,locVars},[],2,'omitnan') ;
    
    % remove those variables from the table
    snp_eqtls_pvalues(:,locVars) = [] ;
    
    % add the variable to the table
    snp_eqtls_pvalues.(char(curVar)) = cur_pvalues ;

end

% get the start of the eqtls agains
eqtl_start = find( ismember(snp_eqtls_pvalues.Properties.VariableNames,...
    'Liver'),true) ;

% remove the rows that have all nans
snp_eqtls_pvalues(all(isnan(snp_eqtls_pvalues{:,eqtl_start:end}),2),:) =[];

clustData = fillmissing(snp_eqtls_pvalues{:,eqtl_start:end},'constant',1);

% now covert the clustData into those which are significant and those with
% are not signficant groups for p-value < 0.05
clustData = double(clustData < 0.05) ;

% produce a clutergram
cgo2 = clustergram(clustData','Colormap',redbluecmap,'RowLabels',...
    snp_eqtls_pvalues.Properties.VariableNames(eqtl_start:end),...
    'ColumnLabels', snp_eqtls_pvalues.Variant ,'Linkage',...
    'complete','ColumnPDist','euclidean','Dendrogram',20);

clear snps_annots tsneData tsneColours biomart cgo cgo2 cgoAfr ...
    cgoAll clustData curVar ans aa b explode ii myTable pieLabel ...
    locVars yValues loc_1kg  gwasAFR gwasEUR eqtl_start ...
    cur_pvalues numOfGenes nonLDsnps nonSigFisher theVars ...
    eurDisease afrDisease finalLabels a eurSigSnps afrSigSnps 

%% Add the snps to the intergated data 

% ***********************************************************************
% ***********************************************************************

% load the sift4g annotations
if ~exist('Intergrated Epigenetic Data.xlsx','file')
    
    % load the snp annotations
    sift4g = snp_annots(:, ...
        {'rs_dbSNP','SIFT4G_pred','SIFT4G_converted_rankscore'} );
    
    % the pie labels have multiple labels get those which have atleast 1
    % predicted pathogenic variant
    sift4g.SIFT4G_pred(contains(sift4g.SIFT4G_pred,'D')) = {'D'} ;
    sift4g.SIFT4G_pred(~contains(sift4g.SIFT4G_pred,'D')) = {'T'} ;
    sift4g.Properties.VariableNames(1) = "Variant" ;
    
    % here are all the datasets
    intData = snpFreq ;
    
    % add the prediction to the table
    [~,locA, locB] = intersect(intData.Variant, sift4g.Variant,'Stable') ;
    intData.sift4(locA) = sift4g.SIFT4G_pred(locB) ;
    intData.sift4_rank(locA) = sift4g.SIFT4G_converted_rankscore(locB) ;
    
    % add the chrom and position to the intData
    epiVariantqc = readtable('variantqc_epi.txt');
    epiVariantqc = unique(epiVariantqc) ;
    epiVariantqc.HugoSymbol = [] ;
    epiVariantqc.Properties.VariableNames(1:2) = {'Chrom','Position'} ;
    intData = innerjoin(intData, epiVariantqc, "LeftKey","Variant", ...
        "RightKey","rsid") ;
    
    % throw in an assertion
    assert(all( ismember( intData.HugoSymbol, intData.nearest_genes)))
    intData.nearest_genes = [] ;
    
    % add the eQTLs information to the table
    [these,locA] = intersect(intData.Variant, ...
        snp_eqtls_pvalues.Variant,'Stable') ;
    intData.eqtl(locA) = true(length(these),1) ;
    
    % add the gwas catalog terms to the integrated dataset
    [~,locA, locB] = intersect(intData.Variant,gwas.Variant,'Stable') ;
    intData.gwasDisease(locA) = gwas.DISEASE_TRAIT(locB) ;
    
    % add the mQTLs infroamtoin to the table and add the variable names that I
    % did not read when I merge the dataset
    mqtl = readtable('epigenetic_mQTL.txt') ;
    mqtl.Properties.VariableNames(11:end) = ...
        {'Status','CpG','beta','t-stat','p-value','mqtl_FDR','SNP_chr',...
        'SNP_site','SNP_alleles','CpG_chr','CpG_site'} ;
    mqtl = mqtl(mqtl.mqtl_FDR < 0.05, :) ;
    [~,locA, locB] = intersect(intData.Variant, mqtl.Variant,'Stable') ;
    intData.mqtl_beta(locA) = mqtl.beta(locB) ;
    intData.mqtl = intData.mqtl_beta ~= 0 ;
    
    % add the hQTL information to the table
    hqtl = readtable('epigenetic_hQTL.txt') ;
    [these,locA, locB] = intersect(intData.Variant, hqtl.Variant,'Stable');
    intData.hqtl_beta(locA) = hqtl.beta(locB) ;
    intData.hqtl_mod(locA) = hqtl.mod(locB) ;
    intData.hqtl(locA) = true(length(these),1) ;
    
    % add the sQTLs to the table
    % add the hQTL information to the table
    sqtl = readtable('epigenetic_sQTL.txt') ;
    sqtl = innerjoin(intData(:,{'Variant','Position','Chrom'}) , ...
        sqtl,"Keys",{'Position','Chrom'} ) ;
    [these,locA, locB] = intersect(intData.Variant, sqtl.Variant,'Stable');
    intData.sqtl_tissue(locA) = sqtl.tissue(locB) ;
    intData.sqtl_fdr(locA) = sqtl.fdr(locB) ;
    intData.sqtl(locA) = true(length(these),1) ;
    
    % remove some unrequired variables
    intData = removevars(intData,{'OddRatio','SigCount'}) ;
    intData = movevars(intData,{'Chrom','Position','ref','alt'},'After',...
        'Variant');
    
    intData = addvars( intData, intData.af_AFR./intData.af_EUR , ...
        'Before','freqDiff','NewVariableNames',{'freqDiv'});
    
    writetable(intData,'Intergrated Epigenetic Data.xlsx')
    
    clear these locA locB hqtl mqtl sqtl sift4g epiVariantqc mutations ...
        ukbFreq  
end

%% Integration with OMMIM, eQTLs, mQTLs, eQTLs, GWAS catalog

% load the intData and and change the highgroup to double 
if ~exist('intData','var')
    intData = readtable('Intergrated Epigenetic Data.xlsx', ...
        'Format','auto');
    intData = intData( ismember(intData.Variant,snpFreq.Variant), :) ;
end
intData.highGroup = categorical(intData.highGroup) ;

% save a dataset for the project member
if ~exist('epigenetics_summary.xlsx','file')
    % make a copy of the data 
    projectData = intData ;

    % remove the unrequired variable 
    projectData = removevars(projectData, ...
        ["sift4","sift4_rank","mqtl_beta","hqtl_beta","hqtl_mod",...
        "sqtl_tissue","sqtl_fdr"]);
    projectData.Properties.VariableNames(14) = "gwasPhenotype";
    projectData = movevars(projectData,"gwasPhenotype","After","FDR") ;

    % save the file to excel 
    writetable(projectData,'epigenetics_summary.xlsx')
    clear projectData
end

% the colors to use the for two groups
groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10  ] ;

% return only the variants that have a disease, sqtl, eqtl, hqtl, and mqtl
plotData = intData( any(intData{:,{'eqtl','mqtl','sqtl','hqtl'}},2), :);

% here are the variable to plot, and scatter shapes, colours
myVars = {'eqtl','mqtl','sqtl','hqtl'} ;
theTitles = {'eQTLs','mQTLs','sQTLs','hQTLs'} ;
myColors = [0 0.447 0.741; 0.494 0.184 0.556;0.635 0.078 0.184; ...
    0.3010 0.7450 0.9330] ;
myShapes = ["o","square","diamond","hexagram"] ;

% here are the letters for plots 
theLetters = 'a':'d';

% here are the figures
figure(200)
tiledlayout(2,2,'TileSpacing','loose'); % ,'padding','compact'

figure(201)
tiledlayout(2,2,'TileSpacing','loose'); % ,'padding','compact'

% plot the figures in a loop 
for ii = 1:length(myVars)

    % plot the graph before removing the causal SNPs
    figure(200)
    nexttile
    hold on

    % plot for the EUR
    scatter(plotData.af_AFR(plotData.highGroup == "EUR"), ...
        plotData.af_EUR(plotData.highGroup == "EUR"), 30,'filled', ...
        'MarkerFaceColor',groupColors(2,:),'Marker','o',...
        'MarkerEdgeColor',groupColors(2,:),'MarkerFaceAlpha',0.3)

    % for the AFR
    scatter(plotData.af_AFR(plotData.highGroup == "AFR"), ...
        plotData.af_EUR(plotData.highGroup == "AFR"), 30,'filled', ...
        'MarkerFaceColor',groupColors(1,:),'Marker','o',...
        'MarkerEdgeColor',groupColors(1,:),'MarkerFaceAlpha',0.3)

    % plot the non significant data
    nonSigFisher = plotData.FDR > 0.05 ;
    scatter(plotData.af_AFR(nonSigFisher), ...
        plotData.af_EUR(nonSigFisher), ...
        30,'filled','MarkerFaceColor',[ 0.5 0.5 0.5 ],'Marker','o',...
        'MarkerEdgeColor',[ 0.5 0.5 0.5 ],'MarkerFaceAlpha',0.3)

    % specific the new colours plot the scatter plots
    scatter(plotData.af_AFR(plotData.(myVars{ii}) == true), ...
        plotData.af_EUR(plotData.(myVars{ii}) == true), 100, ...
        "Marker",myShapes{ii},"MarkerEdgeColor",myColors(ii,:), ...
        "MarkerFaceAlpha",1 ,'LineWidth',2)

    % edit the chart elements
    set(gca,'FontSize',14,'LineWidth',1.5,'Box','off','TickDir','out')
    ylabel('AF in EUR')
    xlabel('AF in AFR');
    legend([{'EUR Sig','AFR sig','Non Sig'},theTitles{ii}], ...
        'Location','best')

    % add the letter to the plot
    text(-0.1, 1 ,theLetters(1),'Units','normalized', ...
        'FontWeight','bold','FontSize',28)

    % hold off the figure
    hold off

    figure(201)
    nexttile
    hold on 
    % specific the new colours plot the scatter plots
    scatter(plotData.af_AFR(plotData.(myVars{ii}) == true), ...
        plotData.af_EUR(plotData.(myVars{ii}) == true), 50, ...
        "Marker",myShapes{ii},"MarkerEdgeColor",myColors(ii,:), ...
        'LineWidth',2)

    % plot the disease variants
    dxData = plotData(plotData.(myVars{ii}) == true,:) ;
    scatter( dxData.af_AFR( ~ismissing(dxData.gwasDisease )), ...
        dxData.af_EUR( ~ismissing(dxData.gwasDisease )) , 50, ...
        "Marker",myShapes{ii},"MarkerEdgeColor",[1 0 0], ...
        "MarkerFaceColor",[1 0 0],'LineWidth',2)

    % here are the number of GWAS phenotypes
    numPheno = sum(~ismissing(dxData.gwasDisease)) ;
    numNoPheno = height(dxData) - numPheno ;

    % add the legend
    lgd = legend({sprintf('%d No',numNoPheno), ...
        sprintf('%d Yes',numPheno)}, 'Location','best');
    title(lgd,'GWAS Phenotype','FontSize',12)

    % edit the chart elements
    set(gca,'FontSize',14,'LineWidth',1.5,'Box','off','TickDir','out')
    ylabel('AF in EUR')
    xlabel('AF in AFR')
    title( [sprintf('%d ',sum( plotData.(myVars{ii}) == true) ), ...
        theTitles{ii}] , 'FontSize',16,'FontWeight','bold')

    % add the letter to the plot
    text(-0.1, 1 ,theLetters(1),'Units','normalized', ...
        'FontWeight','bold','FontSize',28)

    hold off

    % delete the current letter
    theLetters(1) = [];
    
    % print the variant that are significantly more frequent in one group 
    tempData = plotData( plotData.(myVars{ii}) == true, :)  ;
    fprintf('\nThese are the %s variants\n', myVars{ii} )
    summary(tempData.highGroup)
end

clear tempData

%% Plot the distribution of SNP differences 

% here are the group colours again 
groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10  ] ;

% plot the histogram for europeans
figure()

% plot the bar graphs using a tiled layout
tiledlayout(2,1,'padding','compact');
theLetters = ['a','b'] ;

% Plot the data on a bargraph of on the tile
nexttile

% ******************* here is the first histogram *****************
histogram( log10(intData.af_EUR), 100, 'FaceColor', groupColors(2,:) )
set(gca,'FontSize',12,'LineWidth',1,'Box','off','TickDir','out', ...
    'Xlim',[-5.2, 0.2])
hold on 
% plot the histogram of africans
histogram( log10(intData.af_AFR), 100,'FaceColor', groupColors(1,:))

% edit the plots 
xlabel('Log_1_0(SNP Frequency) Bins')
ylabel('Frequency in Bin')
set(gca,'FontSize',12,'LineWidth',1,'Box','off','TickDir','out')
legend({'EUR','AFR'},'Location','best')

% add a letter to the figure and remove the letter from the array
text(-0.07, 1 ,theLetters(1),'Units','normalized','FontWeight','bold', ...
    'FontSize',28)
theLetters(1) = [];

hold off

% save the figure
saveas(gcf,'Bins_of_SNPs.fig','fig')
saveas(gcf,'Bins_of_SNPs.png','png')

% ******************* here is the second histogram *****************
nexttile
hold on 

% plot the histogram of africans
histogram( log10( intData.freqDiv( intData.highGroup == "AFR")) ...
    , 100,'FaceColor', groupColors(1,:))

% plot the histogram for europeans
% divide by one because these values are an inverse of the AFR values 
histogram( log10( 1./intData.freqDiv( intData.highGroup == "EUR")), ...
    100, 'FaceColor', groupColors(2,:) )

% edit the plots 
xlabel('Absolute Log_1_0(AFR AF/ EUR AF) Bins')
ylabel('Frequency in Bin')
set(gca,'FontSize',12,'LineWidth',1,'Box','off','TickDir','out')
legend({'AFR','EUR'},'Location','best')

% add a letter to the figure and remove the letter from the array
text(-0.07, 1 ,theLetters(1),'Units','normalized','FontWeight','bold', ...
    'FontSize',28)
theLetters(1) = [];

hold off

% save the figure
saveas(gcf,'AF_Bins_of_SNPs.fig','fig')
saveas(gcf,'AF_Bins_of_SNPs.png','png')

%% Apply A Fisher Exact Test to find if QTLs are more assoicated with Dx

% preallocate the table and the variable names 
chiTable = table( {'eQTLs';'mQTLs';'sQTLs';'hQTLs';'xQTLs'}, ...
    'VariableNames',{'QTLs'}) ;
myVars = {'eqtl','mqtl','sqtl','hqtl','xqtl'} ;

% here are the actual QTLs to be remove from the table
qtls = {'eqtl','mqtl','sqtl','hqtl'} ;

% run the data in a loop
for ii = 1:height(chiTable)
    
    % get the data for the current qTLS
    if ii < height(chiTable)
        % ge thte current data 
        curData = intData(:,[myVars(ii), 'gwasDisease']) ;

        % remove the rows that have other eQTLs
        curQTLs = qtls;
        curQTLs(ii) = [];
        curData( any( intData{:,curQTLs}, 2) & ...
            intData.(myVars{ii}) ~= true, :) = [] ;

        % plot a word cloud of the traits 
        plotData = curData.gwasDisease( ...
            ~ismissing(curData.gwasDisease) & curData.(myVars{ii}) == true);
        
        % remove levels, count, high and low with ''
        plotData = regexprep(plotData, {'levels','count','high','low', ...
            'volume','total','index','function','mtag', ...
            'mean','time','body','adjusted','type','ratio'},'', ...
            'ignorecase');

        % preprocess the words
        plotData = preprocessText(plotData) ;

        % here the word cloud
        figure()
        wordcloud(plotData)
        title(chiTable.QTLs{ii})

        % now convert to logical array 
        curData.gwasDisease = ~ismissing(curData.gwasDisease) ;

    else
        % get any variants that is an xQTL
        curData = intData(:, {'gwasDisease'}) ;
        curData = addvars(curData, any( ...
            intData{:, {'eqtl','mqtl','sqtl','hqtl'}} , 2),...
            'NewVariableNames',{'xqtl'},'Before',1) ;
        curData.gwasDisease = ~ismissing(curData.gwasDisease) ;
    end


    % get the data for the current fisher exact test 
    fisherData = table('Size',[2,2],'VariableTypes',{'double','double'},...
        'VariableNames',{'Trait','NoTrait'} ) ;
    
    % get the QTL that have gwas phenotypes and those that do not have gwas
    % phenotypes
    fisherData.Trait = [ sum( curData.(myVars{ii}) == true & ...
        curData.gwasDisease == true)  ; ...
        sum( curData.(myVars{ii}) == true & ...
        curData.gwasDisease == false ) ] ;

    % get the non QTL that have gwas phenotypes and those that do not have
    % gwas phenotypes
    fisherData.NoTrait = [ sum( curData.(myVars{ii}) == false & ...
        curData.gwasDisease == true)  ; ...
        sum( curData.(myVars{ii}) == false & ...
        curData.gwasDisease == false ) ] ;

    % now perform the chi-square test
    [~, p, stats ] = fishertest(fisherData);

    % add the number of traits to the table 
    chiTable.QLT_trait(ii) = fisherData.Trait(1)/fisherData.Trait(2) ;
    chiTable.nonQTL_trait(ii) = ...
        fisherData.NoTrait(1)/fisherData.NoTrait(2) ;
    
    % add the results to the table
    chiTable.pValue(ii) = p ; 
    chiTable.OddRatio(ii) = stats.OddsRatio ; 
    chiTable.LowerBound(ii) = stats.ConfidenceInterval(1);
    chiTable.UpperBound(ii) = stats.ConfidenceInterval(2);
    
end  

% add the variable names to the table and then remove the rows with missing
% results and add the adjusted p value to the table 
chiResults = addvars(chiTable, ...
    mafdr(chiTable.pValue,'BHFDR',true),'After','pValue', ...
    'NewVariableNames','FDR') ;
chiResults = sortrows(chiResults,'pValue','ascend');

% save the results to excel 
writetable(chiResults,'Disease QTL Comparision.xlsx')

clear ii curData tbl chi2 chiPvalue labels ci p stats meanAfricans  ...
    meanWhites chiTable ttestTable chiData xValues groups plotName ...
    curCol locFreq africaChiTable

%% VENN DIAGRAM of Intersect QTLs

% here are the ranges in the data 
saveRange = {'A1','B1','C1','D1'} ;

% here are the variable names again excluding the xqtls
myVars = {'eqtl','mqtl','sqtl','hqtl'} ;

% get array for the QTLs
for ii = 1:length(myVars)

    % here is the current data
    curData = intData( intData.(myVars{ii}) == true, {'Variant'}) ;

    % change the variable names
    curData.Properties.VariableNames = theTitles(ii) ;

    % write to an excel file
    writetable(curData,'QTLs_venn_data.xlsx', 'Range',saveRange{ii}) ;

end

clear curData saveRange

%% Get the recomb rate and Enhancer datasets to be used for making plots 

% add the enhancer information to the table 
enhancer = readtable('epigenetic_enhancer.txt') ;

% read the files from the github repository
if ~exist('genetic_map_GRCh37.txt','file')

    % preallocate the genetic map
    geneticMap = [];

    % loop over the chromsomes
    for ii = 1:22

        fprintf('\nGetting genetic map data for chromsome %d\n',ii)

        % get the url where the data are located 
        url = sprintf(['https://github.com/adimitromanolakis/' ...
            'geneticMap-GRCh37/'...
            'raw/master/genetic_map_GRCh37_chr%d.txt.gz'], ii) ;

        % save the file to the computer 
        websave(sprintf('genetic_map_GRCh37_chr%d.txt.gz',ii),url) ;

        % unzip the file 
        gunzip( sprintf('genetic_map_GRCh37_chr%d.txt.gz',ii) ) ;

        % read the file 
        curMap = readtable( sprintf('genetic_map_GRCh37_chr%d.txt',ii) );

        % add column for the chromosome name 
        curMap = addvars( curMap, repmat(ii,height(curMap),1), ...
            'Before',1, 'NewVariableNames', {'Chrom'}) ;

        % add to the growing genetic map table
        geneticMap = [geneticMap ; curMap] ;

        % delete the file from the computer for each chromosome
        system( sprintf('rm -r genetic_map_GRCh37_chr%d.txt', ii) )
    end
    % save the file to the computer 
    writetable(geneticMap,'genetic_map_GRCh37.txt') 
else
    % load the save file
    geneticMap = readtable('genetic_map_GRCh37.txt') ;
end

%% Let See GWAS Phenotypes of QTLs in AFR and EUR 

% I will need to use the manifest to the disease 

% process the LD data
if ~exist('epigenetics_LD_data.txt','file')
   
    % load the LD data
    fprintf('\nLoading the LD data\n')
    ldData = readtable('GWAScat Small PICS2 2021-09-24.txt') ;

    fprintf('\nProcessing the LD data \n')
    % get only the required variables
    ldData = ldData(:,{'IndexSNP','LinkedSNP','x_hg19_','Rsquare'});
    ldData(ldData.Rsquare < 0.4, :) = [] ;

    % get only the data that exist for africans and europeans for the
    % complete datasets and the causal variant dataset just in case
    % something is being missing
    ldData = ldData( ismember(ldData.IndexSNP, intData.Variant) | ...
        ismember(ldData.LinkedSNP, intData.Variant),:) ;

    % merge the two datasets
    ldData = unique(ldData) ;
    ldData.x_hg19_ = [] ;

    % save the results to a file 
    writetable(ldData,'epigenetics_LD_data.txt')
else
    % load the saved file
    ldData = readtable('epigenetics_LD_data.txt') ;
end

% clean up the manifest file 
% *******************************************************************
% CHANGE THE MANIFEST TO PHYSICAL MEASURE TO COMPLETE THE PAPER 
% Please note that p-values are now stored as natural log p-values to
% avoid underflow (i.e., ln P, not -ln P or -log10 P).
% *******************************************************************

if ~exist('clean_manifest.txt','file')

    fprintf('\nCleaning the manifest file\n')

    % read in the table of results location
    manifest = readtable('Pan-UK Biobank phenotype manifest 2.xlsx');
    manifest.trait_type = categorical(manifest.trait_type);
    manifest.description = categorical( manifest.description);
    manifest.phenocode = str2double(manifest.phenocode);

    % return only the data that has african, AFR and EUR population
    manifest = manifest(contains(manifest.pops,{'EUR'}) & ...
        contains(manifest.pops,{'AFR'}),:);

    % also get the data for continuous traits
    % 'Non-cancer illness code, self-reported' and biomarker measures
    manifest = manifest(manifest.trait_type == "biomarkers" | ...
        manifest.description == "Non-cancer illness code, self-reported",:);

    % redo the categories
    manifest.description = categorical( cellstr(manifest.description));
    manifest.trait_type = categorical(cellstr(manifest.trait_type )) ;

    % save the manifest
    writetable(manifest,'clean_manifest.txt')
%     writecell(manifest.filename,'ukb_epi_gwas_filenames.txt')
else
    % load the clean manifest file
    manifest = readtable('clean_manifest.txt') ;
end

% load the variant qc data for epigenetic genes
variant_epi = readtable('variantqc_epi.txt');

% also remove the unnecessary variables
variant_epi = removevars( variant_epi, {'nearest_genes','alt'});

% preallocate the UKB gwas traits 
ukbTraits = [] ;

% define the number of workers on the cluster
% set up a pool of 20 cores
if multiCoreHighMem == true
    try
    % work on the entire datasets 
    parpool('local',numOfCores)
    catch
    end
else
    % *******************  Comment the prototypying code *****************
    manifest = readtable('Pan-UK Biobank phenotype manifest 2.xlsx') ;
    posTry = ismember(manifest.filename, ...
        'continuous-78-both_sexes-irnt.tsv.bgz') ;
    posTry = find(posTry, true) ;
    % ********************************************************************
end

% process the files % change this to a Parallel for Loop
if exist('ukb_traits_epigenetic.txt','file')
    
    % load the processed data
    ukbTraits = readtable('ukb_traits_epigenetic.txt');
else
    % create the data of the file
    for ii = 1:height(manifest)
        
        % get only the high quality variants
        
        % print something to the screen
        fprintf('\n Reading GWAS summary results for %s: # %d of %d\n',...
            manifest.description{ii}, ii, height(manifest) )
        
        % here is the current file names
        curFileName = replace(manifest.filename{ii},'.tsv.bgz','.txt') ;
        
        % try to get the file from AWS
        try
            % check that the file exist with the wrong name
            if exist([manifest.filename{ii},'.txt'],'file')
                movefile([manifest.filename{ii},'.txt'], curFileName )
            end
            
            % check that file exist on the computer
            if ~exist(curFileName,'file')
                
                % load the was from the internet print something to the
                % screen
                fprintf('\nGetting was data %s #%d of %d from AWS\n', ...
                    manifest.description{ii}, ii , height(manifest))
                
                % here is the file name and the location of the data
                filename = manifest.filename{ii} ;
                
                % check if the file has been downloaded
                if ~exist(filename,'file')
                    
                    % here is the file name and the location of the data
                    url = manifest.aws_link{ii};
                    
                    % get the information from AWS
                    websave(filename,url)
                    
                end
                
                % remname the file
                filename_gz = replace( filename,'bgz','gz' );
                movefile(filename, filename_gz)
                
                % unzip the file and untar the file if this is not already
                % done
                fprintf('\nUnzipping the .gz file for %s\n',filename_gz)
                gunzip(filename_gz)
                
                % renem the file
                filename_tsv = extractBefore(filename,'.bgz') ;
                movefile(filename_tsv, curFileName )
                
            end
        catch
            % could not retrive the file
            continue
        end
        
        % ###############################################################
        % USE a datastore to process the results to return only variants of
        % epigentic genes
        
        % check that a smaller file does not exist
        if ~exist(['epi_',curFileName],'file')
            
            % load the gwas results
            fprintf('\n Reading processed gwas results\n')
            
            ds = datastore(curFileName) ;
            
            % Get a number of partitions to parallelize datastore access
            % over the current parallel pool.
            n = numpartitions(ds,gcp);
            
            % time the running of this programme
            tstart = tic ;
            
            % preallocate the all results table
            allResults = [] ;
            
            parfor kk = 1:n
                
                % Partition the datastore and read the data in each part.
                subds = partition(ds,n,kk);
                
                % By calling the read function within a while loop, you can
                % perform intermediate calculations on each subset of data,
                % and then aggregate the intermediate results at the end.
                while hasdata(subds)
                    
                    % read in a smaller dataset: note that this return only
                    % the data that has complete record and those that were
                    % taken at the first instance i.e., have _0_0
                    curResults = read(subds);
                    
                    % chr: Chromosome of the variant. pos: Position of the
                    % variant in GRCh37 coordinates. ref: Reference allele
                    % on the forward strand. alt: Alternate allele (not
                    % necessarily minor allele). Used as effect allele for
                    % GWAS.
                    
                    % change the first variable to chrom so it matches the
                    % variant_epi variable names
                    fprintf(['\n Getting the signficant variants ' ...
                        'for the results\n'])
                    curResults.Properties.VariableNames(1) = "chrom" ;
                    
                    % get only the require common from the gwas results
                    curResults = curResults(:,...
                        {'chrom','pos','ref','alt','beta_AFR', ...
                        'beta_EUR','pval_AFR','pval_EUR'});
                    
                    % get only the significant variants so that the
                    % processing can be faster and keep on original copy of
                    % the file for plot the gwas statatistics using locus
                    % Zoom
                    
                    fprintf('\nMerging GWAS data with variant manifest\n')
                    % merge the datasets
                    curResults = innerjoin(curResults, variant_epi, ...
                        'Keys',{'chrom','pos','ref'} );
                    
                    % ************************************************
                    % convert the p-values to the actual p-values and not
                    % the natural log values
                    curResults.pval_EUR = exp(curResults.pval_EUR);
                    curResults.pval_AFR = exp(curResults.pval_AFR);
                    % ************************************************
                    
                    % have are all the results
                    allResults = [allResults; curResults] ;
                end
            end
            % save the processed file
            writetable(allResults, ['epi_',curFileName])
            
            telapsed = toc(tstart)/60 ;
            fprintf('\n Run time for one file is %d Minutes\n',telapsed)
        end
        
        % ############# end datastore processing #####################
        
        % load the processed file of only epigentic genes
        curResults = readtable(['epi_',curFileName]) ;
        
        % load the current results agains that only have epigeentic genes
        % as curResultsOG
        
        % get only the significant variants in the data
        curResults = curResults( curResults.pval_EUR < 5e-8 | ...
            curResults.pval_AFR < 5e-8, :) ;
        
        % run this in a parrallel for loop because the dataset is large
        % make an array for subseting the data
        ds = datastore('Intergrated Epigenetic Data.xlsx') ;
        
        % Get a number of partitions to parallelize datastore access over
        % the current parallel pool.
        n = numpartitions(ds,gcp);
        
        % time the running of this programme
        tstart = tic ;
        
        % here is the parallel for loop for the subdatasets
        parfor kk = 1:n
            
            % Partition the datastore and read the data in each part.
            subds = partition(ds,n,kk);
            
            % By calling the read function within a while loop, you can
            % perform intermediate calculations on each subset of data, and
            % then aggregate the intermediate results at the end.
            while hasdata(subds)
                
                % read in a smaller dataset: note that this return only the
                % data that has complete record and those that were taken
                % at the first instance i.e., have _0_0
                curIntData = read(subds);
                
                % find the variants in the intergrated data are associated
                % with the traits
                for jj = 1:height(curIntData)
                    
                    % find that variant in the UK biobank results for the
                    % Europeans
                    if ~any( ismember(curResults.rsid, ...
                            curIntData.Variant(jj)))
                        continue
                    else
                        
                        % print something to the screen
                        if rem(jj,100) == 0
                            fprintf(['\nRunning Analysis for SNP %s ' ...
                                'in subdata #%d of %s #%d: number %d of %d\n'],...
                                curIntData.Variant{jj}, kk, ...
                                manifest.description{ii}, ii, ...
                                jj, height(curIntData))
                        end
                        
                        % get teh current row of the data
                        curRow = curIntData(jj,{'HugoSymbol','GeneClass'...
                            'Variant','Chrom','Position','ref','alt', ...
                            'af_AFR','af_EUR','highGroup','FDR', ...
                            'gwasDisease','eqtl','mqtl','hqtl','sqtl'}) ;
                        
                        % get the location of the variant in the ukb gwas
                        % summary statistics dataset
                        locVariant = ismember(curResults.rsid, ...
                            curRow.Variant);
                        
                        % check that the variants does not appear more than
                        % one
                        if sum(locVariant) > 1
                            locVariant = find( locVariant, true, 'first') ;
                        end
                        
                        % here is the snps to be used for replication
                        theSnp = curRow.Variant ;
                        
                        % check the p-values in African and Europeans are
                        % statistically significant
                        if curResults.pval_EUR(locVariant ) < 5e-8 || ...
                                curResults.pval_AFR(locVariant) < 5e-8
                            
                            % add the pvalues to the table
                            curRow.pval_AFR = ...
                                curResults.pval_AFR(locVariant);
                            curRow.pval_EUR  = ...
                                curResults.pval_EUR(locVariant);
                            
                            % also add the trait involved and if its a
                            % biomarker or not -- THIS  WILL COME FROM THE
                            % MENIFEST FILE
                            curRow.UKB_trait_type = ...
                                manifest.trait_type(ii) ;
                            curRow.UKB_description = ...
                                manifest.description(ii) ;
                            
                            % change the trait type to the coding
                            % description
                            if strcmp( ...
                                    manifest.trait_type(ii) ,'categorical')
                                curRow.UKB_description = ....
                                    manifest.coding_description(ii) ;
                            end
                            
                            % check that both p-values are less than 5e-8
                            % for else replicate the resutls in the
                            % population with a non signficant p-value
                            if curResults.pval_EUR(locVariant )> 5e-8 &&...
                                    curResults.pval_AFR(locVariant) < 5e-8
                                
                                % get the summary statistic for europeans
                                % and convert the Pvalue for AFR to P
                                gwasStats = curResults ;
                                gwasStats.P = gwasStats.pval_EUR ;
                                
                                % replicate the results in the AFR group
                                ldsnps = replicateVariant(theSnp, ...
                                    ldData,gwasStats);
                                
                                % add the lead variants to the table
                                if ~isempty(ldsnps)
                                    % only the table is not empty
                                    curRow.rep_pvalue = ldsnps.P(1) ;
                                    curRow.rep_Variant = ldsnps.rsid(1) ;
                                else
                                    % if the table is empty use the default
                                    % p values for AFR and the varinats of
                                    % AFR
                                    curRow.rep_pvalue = curRow.pval_EUR ;
                                    curRow.rep_Variant = curRow.Variant ;
                                end
                                
                                % include a column that show the only EUR
                                % had a significant gwas result
                                curRow.gwasSigGroup = {'AFR'} ;
                                
                            elseif curResults.pval_EUR(locVariant ) ...
                                    < 5e-8 && ...
                                    curResults.pval_AFR(locVariant) > 5e-8
                                % replicate the results in the AFR groups
                                
                                % get the summary statistic for AFR and
                                % convert the Pvalue for AFR to P
                                gwasStats = curResults ;
                                gwasStats.P = gwasStats.pval_AFR ;
                                
                                % replicate the results in the EUR group
                                ldsnps = replicateVariant(theSnp, ...
                                    ldData,gwasStats);
                                
                                % add the lead variants to the table
                                if ~isempty(ldsnps)
                                    % only the table is not empty
                                    curRow.rep_pvalue = ldsnps.P(1) ;
                                    curRow.rep_Variant = ldsnps.rsid(1) ;
                                else
                                    % if the table is empty use the default
                                    % p values for AFR and the varinats of
                                    % AFR
                                    curRow.rep_pvalue = curRow.pval_AFR  ;
                                    curRow.rep_Variant = curRow.Variant ;
                                end
                                
                                % include a column that show the only AFR
                                % had a significant gwas result
                                curRow.gwasSigGroup = {'EUR'} ;
                                
                            else
                                % add the replication variant and the
                                % replication p-value in africans
                                curRow.rep_pvalue = curRow.pval_AFR  ;
                                curRow.rep_Variant = curRow.Variant ;
                                
                                % include a column that show the both EUR
                                % and AFR had significant gwas results
                                curRow.gwasSigGroup = {'Both'} ;
                                
                            end
                            % add the results to the growing table of UKB
                            % gwas traits
                            ukbTraits = [ukbTraits;curRow ] ;
                        end
                    end
                end
            end
        end
    end
    % save the file
    writetable(ukbTraits,'ukb_traits_epigenetic.txt')
end

% get only the variants that are snpFreq 
ukbTraits = ukbTraits( ismember( ukbTraits.Variant, snpFreq.Variant), :) ;

%% Plot the number of AFR and EUR for Each Trait in the Manifest

% here are the data for AFR and EUR only 
manifestPlot = manifest(:, {'description','n_cases_AFR','n_cases_EUR'});
manifestPlot.description = categorical(manifestPlot.description) ;

% plot the bar graph
b = bar( manifestPlot.description, log10(manifestPlot{:,2:end}) );

% edit some plot features 
set(gca,'LineWidth',1,'Box','off','FontSize',12,'TickDir','out')
xlabel('Biomarkers') ;
ylabel('log10 # of cases') ;
title('\bfUKB Cases for Biomarker Traits','FontSize',14)

% add numbers to the top of the plots
for ii = 1:2
    % also change the face color
    b(ii).FaceColor = groupColors(ii,:) ;
end

% add a legend and title to the figure
legend({'AFR','EUR'},'Box','off','Location','northeastoutside')


%% Disease Associated with the eqtls 

% add a new column that has information of the tissue within which a
% partiuclar snp is an eqtl

% get the start of the eqtls agains
eqtl_start = find(ismember(snp_eqtls_pvalues.Properties.VariableNames,...
    'Liver'),true) ;

% here is the new table
snp_eqtls_disease = snp_eqtls_pvalues(:,1:eqtl_start-1) ;

% here are the variables to use
theVars = snp_eqtls_pvalues.Properties.VariableNames(eqtl_start:end) ;

% loop over the vars
for ii = 1:height(snp_eqtls_pvalues) 
    
    % get all the variables that contains that particular snps
    theTissues = theVars( ~isnan(snp_eqtls_pvalues{ii,eqtl_start:end} ) );
    
    % separate the tissue with a semi colom
    if length(theTissues) > 1
       % join the strings
       theTissues = strjoin(strcat(theTissues,';')) ;
       
       % remove the last semi colon
       theTissues = cellstr(theTissues(1:end-1)) ;
    end
    
    % add the variable to the table
    snp_eqtls_disease.Eqtl_Tissues(ii) = theTissues ;

end

% add the gwas catalog data to the eqtl infomraiton
snp_eqtls_disease = innerjoin( snp_eqtls_disease, ...
    unique(gwas(:, {'Variant','DISEASE_TRAIT'} ) ) ) ;

% remove some unrequired variables
snp_eqtls_disease = removevars(snp_eqtls_disease, ...
    {'OddRatio','SigCount'}) ;

%% Internal Functions

% ******************** another internal function ************************
% Helper function to get field value from structure with error checking
function value = get_field(struct, field_name, default_value)
    if isfield(struct, field_name)
        value = struct.(field_name);
    else
        value = default_value;
    end
end

% ******************** another internal function ************************

function documents = preprocessText(textData)

% Erase URLS.
textData = eraseURLs(textData);

% Tokenize.
documents = tokenizedDocument(textData);

% Remove tokens containing digits.
pat = textBoundary + wildcardPattern + digitsPattern + wildcardPattern + textBoundary;
documents = replace(documents,pat,"");

% Convert to lowercase.
documents = lower(documents);

% Remove short words.
documents = removeShortWords(documents,2);

% Remove stop words.
documents = removeStopWords(documents);

end

% **********************************************************************
% ******************* another internal function ************************

% This function trims genes by factor given in the loop It uses three
% functions: genelowvalfilter, genevarfilter and geneentropyfilter

function rnaSeq = filterOutGenes(rnaSeq, numOfGenes)

TotranponseBack = false;
% same times the table is transponsed so reverse that first
try
    genes = rnaSeq.HugoSymbol; % obtain the genes
catch
    % save the original variable name at position 1
    ogVarName1 = rnaSeq.Properties.VariableNames(1) ;
    
    % for the table with genes as variable names
    rnaSeq = transposeTable(rnaSeq);
    rnaSeq.Properties.VariableNames(1) = "HugoSymbol" ;
 
    % set the values of to transponse bback
    TotranponseBack = true;
    
    % get the genes
    genes = rnaSeq.HugoSymbol; % obtain the genes
end
% obtain expression measurements
expression = rnaSeq{:,2:end};
samples = rnaSeq.Properties.VariableNames(2:end) ;


% remove nan values
nanIndices = any(isnan(expression),2);
expression(nanIndices,:) = [];
genes(nanIndices) = [];
numel(genes)

% ======================== Filter Out Genes =====================

% Gene profiling experiments typically include genes that exhibit little
% variation in their profile and are generally not of interest. These genes
% are commonly removed from the data.

% Mask = genevarfilter(Data) calculates the variance for each gene
% expression profile in Data and returns Mask, which identifies the gene
% expression profiles with a variance less than the 10th percentile. Mask
% is a logical vector with one element for each row in Data. The elements
% of Mask corresponding to rows with a variance greater than the threshold
% have a value of 1, and those with a variance less than the threshold are
% 0.

% ****************************** NOTE ********************************* you
% can vary the number of gene returned by editing line 54 ** > 500
% *********************************************************************

% for times = 1:numOfTimes
while length(genes) > numOfGenes
    mask = genevarfilter(expression);
    
    expression = expression(mask,:);
    genes = genes(mask);
    numel(genes)
    
%     % filter out genes below a certain fold change threshold
%     [~,expression,genes] = ...
%         genelowvalfilter(expression,genes,'absval',log2(2));
%     numel(genes)
    
    % filter genes below a certain percentile: VERY POWERFUL discriminant
    [~,expression,genes] = ...
        geneentropyfilter(expression,genes,'prctile',20);
    numel(genes)
end

% finally convert back to a table
rnaSeq = [genes, array2table(expression,'VariableNames', ...
    samples) ] ;
rnaSeq.Properties.VariableNames(1) = "HugoSymbol" ;

% re transponse the table after making the filter
if TotranponseBack  == true
    rnaSeq = transposeTable(rnaSeq);
    rnaSeq.Properties.VariableNames(1) = ogVarName1; 
end


end

%% *************** Some Internal Function *********************
% *************************************************************

function [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% CBREWER - This function produces a colorbrewer table (rgb data) for a
% given type, name and number of colors of the colorbrewer tables. For more
% information on 'colorbrewer', please visit http://colorbrewer2.org/
% 
% The tables were generated from an MS-Excel file provided on the website
% http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html
%
% 
% [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% INPUT:
%   - ctype: type of color table 'seq' (sequential), 'div' (diverging),
%   'qual' (qualitative) - cname: name of colortable. It changes depending
%   on ctype. - ncol:  number of color in the table. It changes according
%   to ctype and
%            cname
%   - interp_method: interpolation method (see interp1.m). Default is
%   "cubic" )
% 
% A note on the number of colors: Based on the original data, there is only
% a certain number of colors available for each type and name of
% colortable. When 'ncol' is larger then the maximum number of colors
% originally given, an interpolation routine is called (interp1) to produce
% the "extended" colormaps.
%
% Example:  To produce a colortable CT of ncol X 3 entries (RGB) of
%           sequential type and named 'Blues' with 8 colors:
%                   CT=cbrewer('seq', 'Blues', 8);
%           To use this colortable as colormap, simply call:
%                   colormap(CT)
% 
%           To see the various colormaps available according to their types
%           and names, simply call: cbrewer()
%
%  This product includes color specifications and designs developed by
%  Cynthia Brewer (http://colorbrewer.org/).
%
% Author: Charles Robert email: tannoudji@hotmail.com Date: 06.12.2011
% ------------------------------ 18.09.2015  Minor fixes, fixed a bug where
% the 'spectral' color table did not appear in the preview

% load colorbrewer data
addpath('/Users/sinkala/Documents/MATLAB/cbrewer/cbrewer/cbrewer')
load(['/Users/sinkala/Documents/MATLAB/cbrewer/cbrewer/cbrewer/' ...
    'colorbrewer.mat'])

% initialise the colormap is there are any problems
colormap=[];
if (~exist('interp_method', 'var'))
    interp_method='cubic';
end

% If no arguments
if (~exist('ctype', 'var') | ~exist('cname', 'var') | ...
        ~exist('ncol', 'var'))
    disp(' ')
    disp('[colormap] = cbrewer(ctype, cname, ncol [, interp_method])')
    disp(' ')
    disp('INPUT:')
    disp(['  - ctype: type of color table *seq* (sequential),...' ...
        ' *div* (divergent), *qual* (qualitative)'])
    disp('  - cname: name of colortable. It changes depending on ctype.')
    disp(['  - ncol:  number of color in the table. ...' ...
        'It changes according to ctype and cname'])
    disp(['  - interp_method:  interpolation method  ...' ...
        '(see interp1.m). Default is "cubic" )'])
    
    disp(' ')
    disp('Sequential tables:')
    z={'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges',...
        'OrRd','PuBu','PuBuGn','PuRd',...
             'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', ...
             'YlOrRd', 'Spectral'};
    disp(z')     
         
    disp('Divergent tables:')
    z={'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'};
    disp(z')
    
    disp(' ')
    disp('Qualitative tables:')
    %getfield(colorbrewer, 'qual')
    z={'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', ...
        'Set2', 'Set3'};
    disp(z')

    plot_brewer_cmap
    return
end

% Verify that the input is appropriate
ctype_names={'div', 'seq', 'qual'};
if (~ismember(ctype,ctype_names))
    disp('ctype must be either: *div*, *seq* or *qual*')
    colormap=[];
    return
end

if (~isfield(colorbrewer.(ctype),cname))
    disp(['The name of the colortable of type *' ctype ...
        '* must be one of the following:'])
    getfield(colorbrewer, ctype)
    colormap=[];
    return
end

if (ncol>length(colorbrewer.(ctype).(cname)))
%     disp(' ')
%     disp('-----------------------------------------------------------')
%     disp(['The maximum number of colors for table *' cname '* is '
%     num2str(length(colorbrewer.(ctype).(cname)))]) disp(['The new
%     colormap will be extrapolated from these '
%     num2str(length(colorbrewer.(ctype).(cname))) ' values'])
%     disp('-----------------------------------------------------------')
%     disp(' ')
    cbrew_init=colorbrewer.(ctype).(cname)...
        {length(colorbrewer.(ctype).(cname))};
    colormap=interpolate_cbrewer(cbrew_init, interp_method, ncol);
    colormap=colormap./255;
    return
end

if (isempty(colorbrewer.(ctype).(cname){ncol}))
    
    while(isempty(colorbrewer.(ctype).(cname){ncol}))
        ncol=ncol+1;
    end        
    disp(' ')
    disp('--------------------------------------------------------------')
    disp(['The minimum number of colors for table *'...
        cname '* is ' num2str(ncol)])
    disp('This minimum value shall be defined as ncol instead')
    disp('--------------------------------------------------------------')
    disp(' ')
end

colormap=(colorbrewer.(ctype).(cname){ncol})./255;

end

%% Another internal function 

function varargout = venn(varargin)

%VENN   Plot 2- or 3- circle area-proportional Venn diagram
%
%  venn(A, I) venn(Z) venn(..., F) venn(..., 'ErrMinMode', MODE) H =
%  venn(...) [H, S] = venn(...) [H, S] = venn(..., 'Plot', 'off') S =
%  venn(..., 'Plot', 'off') [...] = venn(..., P1, V1, P2, V2, ...)
%
%venn(A, I) by itself plots circles with total areas A, and intersection
%area(s) I. For two-circle venn diagrams, A is a two element vector of
%circle areas [c1 c2] and I is a scalar specifying the area of intersection
%between them. For three-circle venn diagrams, A is a three element vector
%[c1 c2 c3], and I is a four element vector [i12 i13 i23 i123], specifiying
%the two-circle intersection areas i12, i13, i23, and the three-circle
%intersection i123.
%
%venn(Z) plots a Venn diagram with zone areas specified by the vector Z.
%For a 2-circle venn diagram, Z is a three element vector [z1 z2 z12] For a
%3-circle venn, Z is a 7 element vector [z1 z2 z3 z12 z13 z23 z123]
%
%venn(..., F) specifies optional optimization options. VENN uses FMINBND to
%locate optimum pair-wise circle distances, and FMINSEARCH to optimize
%overall three-circle alignment. F is a structure with fields specifying
%optimization options for these functions. F may be a two-element array of
%structures, in which case the first structure is used for FMINBND function
%calls, and the second structure is used for FMINSEARCH function calls.
%
%venn(..., 'ErrMinMode', MODE) Used for 3-circle venn diagrams only. MODE
%can be 'TotalError' (default), 'None', or 'ChowRodgers'. When ErrMinMode
%is 'None', the positions and sizes of the three circles are fixed by their
%pairwise-intersections, which means there may be a large amount of error
%in the area of the three- circle intersection. Specifying ErrMinMode as
%'TotalError' attempts to minimize the total error in all four intersection
%zones. The area of the three circles are kept constant in proportion to
%their populations. The 'ChowRodgers' mode uses the the method proposed by
%Chow and Rodgers [Ref. 1] to draw 'nice' three-circle venn diagrams which
%appear more visually representative of the desired areas, although the
%actual areas of the circles are allowed to deviate from requested values.
%
%H = venn(...) returns a two- or three- element vector to the patches
%representing the circles.
%
%[H, S] = venn(...) returns a structure containing descriptive values
%computed for the requested venn diagram. S is a structure with the
%following fields, where C is the number of circles (N = 2 or 3), Z is the
%number of zones (Z = 3 or 7), and I is the number of intersection areas (1
%or 4)
%
% Radius            C-element vector of circle radii
%
% Position          C*2 array of circle centers
%
% ZoneCentroid      Z*2 array of zone centroids (Can be used for labeling)
%
% CirclePop         C-element vector of supplied circle populations.
%                   (I.e., the 'true' circle areas)
%
% CircleArea        C-element of actual circle areas
%
% CircleAreaError   = (CircleArea-CirclePop)/CirclePop
%
% IntersectPop      I-element vector of supplied intersection populations
%                   (I.e., the 'true' intersection areas)
%
% IntersectArea     I-element vector of actual intersection areas
%
% IntersectError    = (IntersectArea-IntersectPop)/IntersectPop
%
% ZonePop           Z-element vector of supplied zone populations. (I.e.
%                   'true' zone areas
%
% ZoneArea          Z-element vector of actual zone areas.
%
% ZoneAreaError     = (ZoneArea-ZonePop)/ZonePop
% 
%
%[H, S] = venn(..., 'Plot', 'off') S = venn(..., 'Plot', 'off')

% Returns a structure of computed values, without plotting the diagram.
% This which can be useful when S is used to draw custom venn diagrams or
% for exporting venn diagram data to another application. When Plot is set
% to off, the handles vector H is returned as an empty array.
% Alternatively, the command

% S = venn(..., 'Plot', 'off) will return only the output structure.

%[...] = venn(..., P1, V1, P2, V2, ...)

%Specifies additional patch settings in standard Matlab parameter/value
%pair syntax. Parameters can be any valid patch parameter. Values for patch
%parameters can either be single values, or a cell array of length
%LENGTH(A), in which case each value in the cell array is applied to the
%corresponding circle in A.

%Examples
%
%   %Plot a simple 2-circle venn diagram with custom patch properties
%   figure, axis equal, axis off A = [300 200]; I = 150;
%   venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black')
%
%   %Compare ErrMinModes A = [350 300 275]; I = [100 80 60 40]; figure
%   subplot(1,3,1), h1 = venn(A,I,'ErrMinMode','None'); axis image,  title
%   ('No 3-Circle Error Minimization') subplot(1,3,2), h2 =
%   venn(A,I,'ErrMinMode','TotalError'); axis image,  title ('Total Error
%   Mode') subplot(1,3,3), h3 = venn(A,I,'ErrMinMode','ChowRodgers'); axis
%   image, title ('Chow-Rodgers Mode') set([h1 h2], 'FaceAlpha', 0.6)
%
%   %Using the same areas as above, display the error optimization at each
%   iteration. Get the output structure. F = struct('Display', 'iter');
%   [H,S] = venn(A,I,F,'ErrMinMode','ChowRodgers','FaceAlpha', 0.6);
%
%   %Now label each zone for i = 1:7
%       text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), ['Zone '
%       num2str(i)])
%   end
%
%See also patch, bar, optimset, fminbdn, fminsearch
%
%Copyright (C) 2008 Darik Gamble, University of Waterloo.
%dgamble@engmail.uwaterloo.ca
%
%References 1. S Chow and P Rodgers. Extended Abstract: Constructing
%Area-Proportional
%   Venn and Euler Diagrams with Three Circles. Presented at Euler Diagrams
%   Workshop 2005. Paris. Available online:
%   http://www.cs.kent.ac.uk/pubs/2005/2354/content.pdf
%
%2. S Chow and F Ruskey. Drawing Area-Proportional Venn and Euler Diagrams.
%   Lecture Notes in Computer Science. 2004. 2912: 466-477.
%   Springer-Verlag. Available online:
%   http://www.springerlink.com/content/rxhtlmqav45gc84q/
%
%3. MP Fewell. Area of Common Overlap of Three Circles. Australian
%Government
%   Department of Defence. Defence Technology and Science Organisation.
%   2006. DSTO-TN-0722. Available online:
%   http://dspace.dsto.defence.gov.au/dspace/bitstream/1947/4551/4/...
% DSTO-TN-0722.PR.pdf


%Variable overview
%   A0, A   Desired and actual circle areas
%               A = [A1 A2] or [A1 A2 A3]
%   I0, I   Desired and actual intersection areas
%               I = I12 or [I12 I13 I23 I123]
%   Z0, Z   Desired and actual zone areas
%               Z = [Z1 Z2 Z12] or [Z1 Z2 Z3 Z12 Z13 Z23 Z123]
%   x, y    Circle centers
%               x = [x1 x2] or [x1 x2 x3]
%   r       Circle radii
%               r = [r1 r2] or [r1 r2 r3]
%   d       Pair-wise distances between circles
%               d = d12 or [d12 d13 d23]
    


    %Parse input arguments and preallocate settings
    [A0, I0, Z0, nCirc, fminOpts, vennOpts, patchOpts] = ...
        parseArgsIn(varargin);
    [d, x, y, A, I, Z] = preallocVectors (nCirc);
    zoneCentroids = []; %Will only be calculated if needed
    
    %Circle Radii
    r = sqrt(A0/pi);

    %Determine distance between first circle pair
    d(1) = circPairDist(r(1), r(2), I0(1), fminOpts(1));
    
    %Position of second circle is now known
    x(2) = d(1); 
    
    %First intersection area
    I(1) = areaIntersect2Circ(r(1), r(2), d(1));
    
    if nCirc==3
        %Pairwise distances for remaining pairs 1&3 and 2&3
        d(2) = circPairDist(r(1), r(3), I0(2), fminOpts(1)); %d13
        d(3) = circPairDist(r(2), r(3), I0(3), fminOpts(1)); %d23

        %Check triangle inequality
        srtD = sort(d);
        if ~(srtD(end)<(srtD(1)+srtD(2)))
            error('venn:triangleInequality', ...
                'Triangle inequality not satisfied')
        end

        %Guess the initial position of the third circle using the law of
        %cosines
        alpha = acos( (d(1)^2 + d(2)^2 - d(3)^2)  / (2 * d(1) * d(2)) );
        x(3) = d(2)*cos(alpha);
        y(3) = d(2)*sin(alpha);

        %Each pair-wise intersection fixes the distance between each pair
        %of circles, so technically there are no degrees of freedom left in
        %which to adjust the three-circle intersection. We can either try
        %moving the third circle around to minimize the total error, or
        %apply Chow-Rodgers
        
        switch vennOpts.ErrMinMode
            case 'TotalError'
                %Minimize total intersection area error by moving the third
                %circle
                pos = fminsearch(@threeCircleAreaError, [x(3) y(3)],...
                    fminOpts(2));
                x(3) = pos(1);
                y(3) = pos(2);
            case 'ChowRodgers'
                %note that doChowRodgersSearch updates x and y in this
                %workspace as a nested fcn
                doChowRodgersSearch;
        end

        %Make sure everything is 'up to date' after optimization
        update3CircleData;
        
    end
    
    % Are we supposed to plot?
    if vennOpts.Plot
        if isempty(vennOpts.Parent)
            vennOpts.Parent = gca;
        end
        hVenn = drawCircles(vennOpts.Parent, x, y, r, ...
            patchOpts.Parameters, patchOpts.Values);
    else
        hVenn = [];
    end
    
    %Only determine zone centroids if they're needed Needed for output
    %structure
    nOut = nargout;
    if (nOut==1 && ~vennOpts.Plot) || nOut==2
        if nCirc == 2
            %Need to calculate new areas
            A = A0; %Areas never change for 2-circle venn
            Z = calcZoneAreas(2, A, I);
            zoneCentroids = zoneCentroids2(d, r, Z);
        else
            zoneCentroids = zoneCentroids3(x, y, d, r, Z);
        end
    end
        
    %Figure out output arguments
    if nOut==1
        if vennOpts.Plot
            varargout{1} = hVenn;
        else
            varargout{1} = getOutputStruct;
        end
    elseif nOut==2
        varargout{1} = hVenn;
        varargout{2} = getOutputStruct;
    end
        
    
    
    function err = threeCircleAreaError (pos)
        
        x3 = pos(1);
        y3 = pos(2);
        
        %Calculate distances
        d(2) = sqrt(x3^2 + y3^2); %d13
        d(3) = sqrt((x3-d(1))^2 + y3^2); %d23
        
        %Calculate intersections Note: we're only moving the third circle,
        %so I12 is not changing
        I(2:3) = areaIntersect2Circ (r(1:2), r([3 3]), d(2:3)); %I13 andI23
        I(4) = areaIntersect3Circ (r, d); %I123
        
        %Replace 0 (no intersection) with infinite error
        I(I==0) = Inf;
        
        %Error
        err = sum(abs((I-I0)./I0));
        
    end


    function doChowRodgersSearch
        
        %Adapted from Ref. [1]
        
        %Initialize an index matrix to select all 7choose2 zone pairs (21
        %pairs)
        idx = nchoosek(1:7, 2);
                
        %Which zone-zone pairs are considered equal? Zones within 10% of
        %each other considered equal
        zonePairAreas0 = Z0(idx);
        
        %Percent difference in population between the two members of a pair
        ar0 = 2*abs(zonePairAreas0(:,1)...
            -zonePairAreas0(:,2))./sum(zonePairAreas0, 2)*100;
        eqPairCutoff = 10;  
        pairIsEq = ar0<=eqPairCutoff;
        
        %Calculate allowable range for pairs of zones considered unequal
        if any(~pairIsEq)
            %Sort zone areas
            [zUneqAreas0, zUneqAreasSrtIdx] = ...
                sort(zonePairAreas0(~pairIsEq,:), 2);
            
            %Make a real index array out of the inconvenient index sort
            %returns
            n = sum(~pairIsEq);
            zUneqAreasSrtIdx = sub2ind([n,2],[1:n; 1:n]',zUneqAreasSrtIdx);
            
            %rp = (largepopulation/smallpopulation)-1
            rp = zUneqAreas0(:,2)./zUneqAreas0(:,1)-1;
            rpMin = 1 + 0.3*rp;
            rpMax = 1 + 2*rp;
        end
        
        %Preallocate zone error vector
        zoneErr = zeros(1,21); 

        %Initialize independent parameters to search over
        guessParams = [r(1) x(2) r(2) x(3) y(3) r(3)];
        
        %Search!
        pp = fminsearch(@chowRodgersErr, guessParams, fminOpts(2));  
        
        [r(1) x(2) r(2) x(3) y(3) r(3)] = deal(pp(1), pp(2), pp(3),...
            pp(4), pp(5), pp(6));
        
        
        function err = chowRodgersErr (p)
            
            %params = [x2 r2 x3 y3 r3]
            [r(1), x(2), r(2), x(3), y(3), r(3)] = deal(p(1), p(2),...
                p(3), p(4), p(5), p(6));
             
            %After changing x2, r2, x3, y3, and r3, update circle areas,
            %distances, intersection areas, zone areas
            update3CircleData;

            if any(pairIsEq)
                %For zone pairs considered equal, error is equal to square
                %of the distance beyond the cutoff; 0 within cutoff
                zAreas = Z(idx(pairIsEq,:));
                ar = 2*abs(zAreas(:,1)-zAreas(:,2))./sum(zAreas, 2)*100;
                isWithinRange = ar<eqPairCutoff;
                ar(isWithinRange) = 0;
                ar(~isWithinRange) = ar(~isWithinRange) - eqPairCutoff;

                %Amplify error for equal zones with unequal areas
                eqZoneUneqAreaErrorGain = 10;
                ar(~isWithinRange) = ar(~isWithinRange)*...
                    eqZoneUneqAreaErrorGain;

                zoneErr(pairIsEq) = ar.^2;
            end

            if any(~pairIsEq)
                %For zone pairs considered unequal, error is equal to
                %square of the distance from the allowable range of rp

                %rp = (largepopulation/smallpopulation)-1
                zUneqPairAreas = Z(idx(~pairIsEq,:));
                
                %Sort based on the population sizes (determined by parent
                %function doChowRodgersSearch)
                zUneqPairAreas = zUneqPairAreas(zUneqAreasSrtIdx);
                rp = zUneqPairAreas(:,2)./zUneqPairAreas(:,1)-1;

                lessThanMin = rp<rpMin;
                moreThanMax = rp>rpMax;
                rp(~lessThanMin & ~moreThanMax) = 0;
                
                %Determine how far out of range errors are
                rp(lessThanMin) = rp(lessThanMin) - rpMin(lessThanMin);
                rp(moreThanMax) = rp(moreThanMax) - rpMax(moreThanMax);          

                %Consider the case where rp < rpMin to be more erroneous
                %than the case where rp > rpMax
                tooSmallErrorGain = 10;
                rp(lessThanMin) = rp(lessThanMin)*tooSmallErrorGain;

                zoneErr(~pairIsEq) = rp.^2;
            end
            
            %Total error
            err = sum(zoneErr);
            
        end %chowRodgersErr
        
    end %doChowRodgersSearch

    function update3CircleData
        
        %Circle areas
        A = pi*r.^2;

        %Calculate distances
        d(1) = abs(x(2)); %d12
        d(2) = sqrt(x(3)^2 + y(3)^2); %d13
        d(3) = sqrt((x(3)-d(1))^2 + y(3)^2); %d23

        %Calculate actual intersection areas
        I(1:3) = areaIntersect2Circ (r([1 1 2]), r([2 3 3]), d);%I12I13,I23
        I(4) = areaIntersect3Circ (r, d); %I123

        %Calculate actual zone areas
        Z = calcZoneAreas(3, A, I);

    end

    function S = getOutputStruct
               
        S = struct(...
            'Radius'                ,r                      ,...        
            'Position'              ,[x' y']                ,...
            'ZoneCentroid'          ,zoneCentroids          ,...
            'CirclePop'             ,A0                     ,...
            'CircleArea'            ,A                      ,...
            'CircleAreaError'       ,(A-A0)./A0             ,...
            'IntersectPop'          ,I0                     ,...
            'IntersectArea'         ,I                      ,...
            'IntersectError'        ,(I-I0)./I0             ,...
            'ZonePop'               ,Z0                     ,...
            'ZoneArea'              ,Z                      ,...
            'ZoneAreaError'         ,(Z-Z0)./Z0             );  
        end

end %venn

        
function D = circPairDist (rA, rB, I, opts)
    %Returns an estimate of the distance between two circles with radii rA
    %and rB with area of intersection I opts is a structure of FMINBND
    %search options
    D = fminbnd(@areadiff, 0, rA+rB, opts);
    function dA = areadiff (d)
        intersectArea = areaIntersect2Circ (rA, rB, d);
        dA = abs(I-intersectArea)/I;
    end
end

function hCirc = drawCircles(hParent, xc, yc, r, P, V,c)

    hAx = ancestor(hParent, 'axes');
    nextplot = get(hAx, 'NextPlot');
    
    %P and V are cell arrays of patch parameter/values
    xc = xc(:); yc = yc(:);     %Circle centers
    r = r(:);                   %Radii
    n = length(r);              
    
    %Independent parameter
    dt = 0.05;
    t = 0:dt:2*pi;

    % Origin centered circle coordinates
    X = r*cos(t);
    Y = r*sin(t);
    
    hCirc = zeros(1,n);
    % c = [0.00,0.45,0.74; 0.85,0.33,0.10; 0.47,0.67,0.19]; c =
    % [0,0.450,0.740;0.640,0.080,0.18] ; %blue and orange-red
    c = [ 0.85,0.33,0.10; 0.47,0.67,0.19  ]; % green and orange-brown
    
    % c = {'r', 'g', 'b'};                        %default colors
    fa = {0.6, 0.6, 0.6};                       %default face alpha
    tag = {'Circle1', 'Circle2', 'Circle3'}; 	%default tag
    
    for i = 1:n
        xx = X(i,:)+xc(i);  
        yy = Y(i,:)+yc(i);
        
        if iscell(c)
            % if the color is cell array
            hCirc(i) = patch (xx, yy, c{i}, 'FaceAlpha', fa{i},...
                'Parent', hParent, 'Tag', tag{i});
        else % or my colors
            hCirc(i) = patch (xx, yy, c(i,:), 'FaceAlpha', fa{i},...
                'Parent', hParent, 'Tag', tag{i} ,'LineWidth',1);
        end
        
        if i==1
            set(hAx, 'NextPlot', 'add');
        end
    end
    set(hAx, 'NextPlot', nextplot);
    
    % make the axis invisible
    set(findobj(gcf, 'type','axes'), 'Visible','off')

%     % Custom patch parameter values if ~isempty(P)
% 
%         c = cellfun(@iscell, V);
% 
%         %Scalar parameter values -- apply to all circles if any(~c)
%             set(hCirc, {P{~c}}, {V{~c}});
%         end
% 
%         %Parameters values with one value per circle if any(c)
%             %Make sure all vals are column cell arrays V = cellfun(@(val)
%             (val(:)), V(c), 'UniformOutput', false); set(hCirc, {P{c}},
%             [V{:}])
%         end
%     end
    
end %plotCircles

 
function A = areaIntersect2Circ (r1, r2, d)
    %Area of Intersection of 2 Circles Taken from [2]
    
    alpha = 2*acos( (d.^2 + r1.^2 - r2.^2)./(2*r1.*d) );
    beta  = 2*acos( (d.^2 + r2.^2 - r1.^2)./(2*r2.*d) );
    
    A =    0.5 * r1.^2 .* (alpha - sin(alpha)) ...  
         + 0.5 * r2.^2 .* (beta - sin(beta));
    
end

function [A, x, y, c, trngArea] = areaIntersect3Circ (r, d)
    %Area of common intersection of three circles This algorithm is taken
    %from [3].
    %   Symbol    Meaning
    %     T         theta p         prime pp        double prime
        
    %[r1 r2 r3] = deal(r(1), r(2), r(3)); [d12 d13 d23] = deal(d(1), d(2),
    %d(3));

    %Intersection points
    [x,y,sinTp,cosTp] = intersect3C (r,d);
    
    if any(isnan(x)) || any(isnan(y))
        A = 0;
        %No three circle intersection
        return
    end
    
    %Step 6. Use the coordinates of the intersection points to calculate
    %the chord lengths c1, c2, c3:
    i1 = [1 1 2];
    i2 = [2 3 3];
    c = sqrt((x(i1)-x(i2)).^2 + (y(i1)-y(i2)).^2)';

    %Step 7: Check whether more than half of circle 3 is included in the
    %circular triangle, so as to choose the correct expression for the area
    lhs = d(2) * sinTp;
    rhs = y(2) + (y(3) - y(2))/(x(3) - x(2))*(d(2)*cosTp - x(2));
    if lhs < rhs
        sign = [-1 -1 1];
    else
        sign = [-1 -1 -1];
    end
    
    %Calculate the area of the three circular segments.
    ca = r.^2.*asin(c/2./r) + sign.*c/4.*sqrt(4*r.^2 - c.^2);

    trngArea = 1/4 * sqrt( (c(1)+c(2)+c(3))*(c(2)+c(3)-c(1))*...
        (c(1)+c(3)-c(2))*(c(1)+c(2)-c(3)) );
    A = trngArea + sum(ca);
    
end

function [x, y, sinTp, cosTp] = intersect3C (r, d)
    %Calculate the points of intersection of three circles Adapted from
    %Ref. [3]
    
    %d = [d12 d13 d23] x = [x12; x13; x23] y = [y12; y13; y23]

    %   Symbol    Meaning
    %     T         theta p         prime pp        double prime
    
    x = zeros(3,1);
    y = zeros(3,1);
     
    %Step 1. Check whether circles 1 and 2 intersect by testing d(1)
    if ~( ((r(1)-r(2))<d(1)) && (d(1)<(r(1)+r(2))) )
        %x = NaN; y = NaN; bigfix: no returned values for sinTp, cosTp
        [x, y, sinTp, cosTp] = deal(NaN);
        return
    end

    %Step 2. Calculate the coordinates of the relevant intersection point
    %of circles 1 and 2:
    x(1) = (r(1)^2 - r(2)^2 + d(1)^2)/(2*d(1));
    y(1) = 0.5/d(1) * sqrt( 2*d(1)^2*(r(1)^2 + r(2)^2) - (r(1)^2 ...
        - r(2)^2)^2 - d(1)^4 );

    %Step 3. Calculate the values of the sines and cosines of the angles tp
    %and tpp:
    cosTp  =  (d(1)^2 + d(2)^2 - d(3)^2) / (2 * d(1) * d(2));
    cosTpp = -(d(1)^2 + d(3)^2 - d(2)^2) / (2 * d(1) * d(3));
    sinTp  =  (sqrt(1 - cosTp^2));
    sinTpp =  (sqrt(1 - cosTpp^2));

    %Step 4. Check that circle 3 is placed so as to form a circular
    %triangle.
    cond1 = (x(1) - d(2)*cosTp)^2 + (y(1) - d(2)*sinTp)^2 < r(3)^2;
    cond2 = (x(1) - d(2)*cosTp)^2 + (y(1) + d(2)*sinTp)^2 > r(3)^2;
    if  ~(cond1 && cond2)
        x = NaN; y = NaN;
        return
    end

    %Step 5: Calculate the values of the coordinates of the relevant
    %intersection points involving circle 3
    xp13  =  (r(1)^2 - r(3)^2 + d(2)^2) / (2 * d(2));
    %yp13  = -0.5 / d(2) * sqrt( 2 * d(2)^2 * (r(2)^2 + r(3)^2) - (r(1)^2 -
    %r(3)^2)^2 - d(2)^4 );
    yp13  = -0.5 / d(2) * sqrt( 2 * d(2)^2 * (r(1)^2 + r(3)^2) ...
        - (r(1)^2 - r(3)^2)^2 - d(2)^4 );

    x(2)   =  xp13*cosTp - yp13*sinTp;
    y(2)   =  xp13*sinTp + yp13*cosTp;

    xpp23 =  (r(2)^2 - r(3)^2 + d(3)^2) / (2 * d(3));
    ypp23 =  0.5 / d(3) * sqrt( 2 * d(3)^2 * (r(2)^2 + r(3)^2) ...
        - (r(2)^2 - r(3)^2)^2 - d(3)^4 );

    x(3) = xpp23*cosTpp - ypp23*sinTpp + d(1);
    y(3) = xpp23*sinTpp + ypp23*cosTpp;

end



function z = calcZoneAreas(nCircles, a, i)
    
    %Uses simple set addition and subtraction to calculate the zone areas
    %with circle areas a and intersection areas i

    if nCircles==2
        %a = [A1 A2] i = I12 z = [A1-I12, A2-I12, I12]
        z = [a(1)-i, a(2)-i, i];
    elseif nCircles==3
        %a = [A1  A2  A3] i = [I12 I13 I23 I123] z = [A1-I12-I13+I123,
        %A2-I12-I23+I123, A3-I13-I23+I123, ...
        %     I12-I123, I13-I123, I23-I123, I123];
        z = [a(1)-i(1)-i(2)+i(4), a(2)-i(1)-i(3)+i(4), ...
            a(3)-i(2)-i(3)+i(4), ...
                i(1)-i(4), i(2)-i(4), i(3)-i(4), i(4)];
    else
        error('')
        %This error gets caught earlier in the stack w. better error msgs
    end
end

function [Cx, Cy, aiz] = centroid2CI (x, y, r)

    %Finds the centroid of the area of intersection of two circles.
    %Vectorized to find centroids for multiple circle pairs x, y, and r are
    %nCirclePairs*2 arrays Cx and Cy are nCirclePairs*1 vectors

    %Centroid of the area of intersection of two circles
    n = size(x,1);
    xic = zeros(n,2);
    az = zeros(n,2);
    
    dx = x(:,2)-x(:,1);
    dy = y(:,2)-y(:,1);
    d = sqrt(dx.^2 + dy.^2);
    
    %Translate the circles so the first is at (0,0) and the second is at
    %(0,d) By symmetry, all centroids are located on the x-axis. The two
    %circles intersect at (xp, yp) and (xp, -yp)
    xp = 0.5*(r(:,1).^2 - r(:,2).^2 + d.^2)./d;

    %Split the inner zone in two Right side (Area enclosed by circle 1 and
    %the line (xp,yp) (xp,-yp) Angle (xp,yp) (X1,Y1) (xp,-yp)
    alpha = 2*acos(xp./r(:,1));
    %Area and centroid of the right side of the inner zone
    [xic(:,1) az(:,1)] = circleChordVals (r(:,1), alpha);
    %Angle (xp,yp) (X2,Y2) (xp,-yp)
    alpha = 2*acos((d-xp)./r(:,2));
    %Area and centroid of the left side of the inner zone
    [xic(:,2) az(:,2)] = circleChordVals (r(:,2), alpha);
    xic(:,2) = d - xic(:,2);

    %Thus the overall centroid  & area of the inner zone
    aiz = sum(az,2);
    Cx = sum(az.*xic,2)./aiz;
    
    %Now translate the centroid back based on the original positions of the
    %circles
    theta = atan2(dy, dx);
    Cy = Cx.*sin(theta) + y(:,1);
    Cx = Cx.*cos(theta) + x(:,1);
    
end

function centroidPos = zoneCentroids2 (d, r, Z)
    
    centroidPos = zeros(3,2);
    
    %Find the centroids of the three zones in a 2-circle venn diagram By
    %symmetry, all centroids are located on the x-axis. First, find the
    %x-location of the middle (intersection) zone centroid
    
    %Centroid of the inner zone
    centroidPos(3,1) = centroid2CI([0 d], [0 0], r);
    
    %Now, the centroid of the left-most zone is equal to the centroid of
    %the first circle (0,0) minus the centroid of the inner zone
    centroidPos(1,1) = -centroidPos(3,1)*Z(3)/Z(1);
    
    %Similarly for the right-most zone; the second circle has centroid at
    %x=d
    centroidPos(2,1) = (d*(Z(2)+Z(3)) - centroidPos(3,1)*Z(3))/Z(2);
    
end

function centroidPos = zoneCentroids3 (x0, y0, d, r, Z)

    Z = Z(:);
        
    %Get area, points of intersection, and chord lengths
    [act, xi, yi, c, atr] = areaIntersect3Circ (r, d);
    atr = atr(:);
    r = r(:);
    
    %Area and centroid of the triangle within the circular triangle is
    xtr = sum(xi/3); 
    ytr = sum(yi/3);
        
    %Now find the centroids of the three segments surrounding the triangle
    i = [1 2; 1 3; 2 3]; 
    xi = xi(i); yi = yi(i);
    [xcs, ycs, acs] = circSegProps (r(:), x0(:), y0(:), xi, yi, c(:));
    
    %Overall centroid of the circular triangle
    xct = (xtr*atr + sum(xcs.*acs))/act;
    yct = (ytr*atr + sum(ycs.*acs))/act;
    
    %Now calculate the centroids of the three two-pair intersection zones
    %(Zones 12 13 23) Entire zone centroid/areas

    %x, y, and r are nCirclePairs*2 arrays Cx and Cy are nCirclePairs*1
    %vectors
    i = [1 2; 1 3; 2 3];
    [x2c, y2c, a2c] = centroid2CI (x0(i), y0(i), r(i));
    
    %Minus the three-circle intersection zone
    xZI2C = (x2c.*a2c - xct*act)./(a2c-act);
    yZI2C = (y2c.*a2c - yct*act)./(a2c-act);
    
    x0 = x0(:);
    y0 = y0(:);
    
    %Finally, the centroids of the three circles minus the intersection
    %areas
    i1 = [4 4 5]; i2 = [5 6 6];
    j1 = [1 1 2]; j2 = [2 3 3];
    x1C = (x0*pi.*r.^2 - xZI2C(j1).*Z(i1) - xZI2C(j2).*Z(i2) - ...
        xct*act)./Z(1:3);
    y1C = (y0*pi.*r.^2 - yZI2C(j1).*Z(i1) - yZI2C(j2).*Z(i2) - ...
        yct*act)./Z(1:3);
    
    %Combine and return
    centroidPos = [x1C y1C; xZI2C yZI2C; xct yct];
end


function [x, a] = circleChordVals (r, alpha)
    %For a circle centered at (0,0), with angle alpha from the x-axis to
    %the intersection of the circle to a vertical chord, find the
    %x-centroid and area of the region enclosed between the chord and the
    %edge of the circle adapted from
    %http://mathworld.wolfram.com/CircularSegment.html
    a = r.^2/2.*(alpha-sin(alpha));                         %Area
    x = 4.*r/3 .* sin(alpha/2).^3 ./ (alpha-sin(alpha));    %Centroid
end

function [xc, yc, area] = circSegProps (r, x0, y0, x, y, c)

    %Translate circle to (0,0)
    x = x-[x0 x0];
    y = y-[y0 y0];

    %Angle subtended by chord
    alpha = 2*asin(0.5*c./r);
       
    %adapted from http://mathworld.wolfram.com/CircularSegment.html
    area = r.^2/2.*(alpha-sin(alpha));                         %Area
    d   = 4.*r/3 .* sin(alpha/2).^3 ./ (alpha-sin(alpha));    %Centroid
   
    %Perpindicular bisector of the chord
    m = -(x(:,2)-x(:,1))./(y(:,2)-y(:,1));
    
    %angle of bisector
    theta = atan(m);
    
    %centroids
    xc = d.*cos(theta);
    yc = d.*sin(theta);
    
    %Make sure we're on the correct side Point of intersection of the perp.
    %bisector and the circle perimeter
    xb = (x(:,1)+x(:,2))/2;
    xc(xb<0) = xc(xb<0)*-1;
    yc(xb<0) = yc(xb<0)*-1;
    
    %Translate back
    xc = xc + x0;
    yc = yc + y0;
end


function [A0, I0, Z0, nCircles, fminOpts, vennOpts, patchOpts] = ...
    parseArgsIn (args)

    [A0, I0, Z0] = deal([]);
    nIn = length(args);
    badArgs = false;
    
    %Get the easy cases out of the way
    if nIn == 0
        badArgs = true;
    elseif nIn == 1
        %venn(Z)
        Z0 = args{1};
        nIn = 0;
    elseif nIn == 2
        if isnumeric(args{2})
            %venn (A,I)
            [A0, I0] = deal(args{1:2});
            nIn = 0;
        else
            %venn (Z, F)
            Z0 = args{1};
            args = args(2);
            nIn = 1;
        end
    else
        %Find the first non-numeric input arg
        i = find(~cellfun(@isnumeric, args), 1);
        if i == 2
            %venn(Z, ....)
            Z0 = args{1};
        elseif i == 3
            %venn(A, I, ...)
            [A0, I0] = deal(args{1:2});
        else
            badArgs = true;
        end
        nIn = nIn - i + 1;
        args = args(i:end);
    end
    
    if badArgs
        error('venn:parseInputArgs:unrecognizedSyntax',...
            'Unrecogized input syntax')
    end
    try
        [A0, I0, Z0] = parseInputAreas (A0, I0, Z0);
    catch
        error('venn:parseArgsIn:parseInputAreas',...
            'Incorrect size(s) for area vector(s)')
    end
    nCircles = length(A0);
    nZones = length(Z0);
             
    %Any arguments left?
    if nIn > 0 
        
        if isstruct(args{1})
            %FMIN search options
            f = args{1};
           
            nIn = nIn - 1;
            if nIn>0, args = args(2:end); end

            if length(f) == 1
                %Just double up
                fminOpts = [f f];
            elseif length(f) == 2
                %ok
                fminOpts = f;
            else
                error('venn:parseArgsIn',...
                    'FMINOPTS must be a 1 or 2 element structure array.')
            end
        else
            %Use defaults
            fminOpts = [optimset('fminbnd'), optimset('fminsearch')];
        end
    else
        %Use defaults
        fminOpts = [optimset('fminbnd'), optimset('fminsearch')];
    end

    %If there's an even number of args in remaining
    if nIn>0 
        if mod(nIn, 2)==0
            %Parameter/Value pairs
            p = args(1:2:end);
            v = args(2:2:end);
            [vennOpts, patchOpts] = parsePVPairs (p, v, nZones);
        else
            error('venn:parseArgsIn', ...
                'Parameter/Value options must come in pairs')
        end
    else
        vennOpts = defaultVennOptions;
        patchOpts = struct('Parameters', [], 'Values', []);
    end

end %parseArgsIn

function [vennOpts, patchOpts] = parsePVPairs (p, v, nZones)

    p = lower(p);

    %Break up P/V list into Venn parameters and patch parameters
    vennParamNames = {'plot', 'errminmode', 'parent'};
    [isVennParam, idx] = ismember(p, vennParamNames);
    idx = idx(isVennParam);
    %vennParams = p(isVennParam);
    vennVals = v(isVennParam);
    
    %First do Patch options
    patchOpts.Parameters = p(~isVennParam);
    patchOpts.Values = v(~isVennParam);
        
    %Now do Venn options
    vennOpts = defaultVennOptions;
        
    %PLOT
    i = find(idx==1, 1);
    if i
        plot = lower(vennVals{i});
        if islogical(plot)
            vennOpts.Plot = plot;
        else
            if ischar(plot) && any(strcmp(plot, {'on', 'off'}))
                vennOpts.Plot = strcmp(plot, 'on');
            else
                error('venn:parsePVPairs', ['Plot must be ''on'',...' ...
                    ' ''off'', or a logical value.'])
            end
        end
    end
    
    %ERRMINMODE
    i = find(idx==2, 1);
    if i
        mode = lower(vennVals{i});
        okModes = {'None', 'TotalError', 'ChowRodgers'};
        [isOkMode, modeIdx] = ismember(mode, lower(okModes));
        if isOkMode                
            vennOpts.ErrMinMode = okModes{modeIdx};
        else
            error('venn:parsePVPairs', ['ErrMinMode must be None,...' ...
                ' TotalError, or ChowRodgers'])
        end
    end

    %PARENT
    i = find(idx==5, 1);
    if i
        h = v{i};
        if length(h)==1 && ishandle(h) 
            vennOpts.Parent = h;
        else
            error('venn:parsePVPairs', ...
                'Parent must be a valid scalar handle')
        end
    end
    
       
end %parsePVPairs

function [A0, I0, Z0] = parseInputAreas (A0, I0, Z0)

    %Switch to row vectors
    A0 = A0(:)';
    I0 = I0(:)';
    Z0 = Z0(:)';

    if isempty(Z0)
        %A0 and I0 supplied
        
        Z0 = calcZoneAreas (length(A0), A0, I0);
    else
        %Z0 supplied
        switch length(Z0)
            case 3
                A0 = Z0(1:2)+Z0(3);
                I0 = Z0(3);
            case 7
                A0 = Z0(1:3)+Z0([4 4 5])+Z0([5 6 6])+Z0(7);
                I0 = [Z0(4:6)+Z0(7) Z0(7)];
            otherwise
                error('')
        end
    end
end

function vennOpts = defaultVennOptions 
    
    vennOpts = struct(...
        'Plot'          ,true               ,...
        'Labels'        ,[]                 ,...
        'PopLabels'     ,false              ,...
        'DrawLabels'    ,false              ,...
        'Parent'        ,[]                 ,...
        'Offset'        ,[0 0]              ,...
        'ErrMinMode'    ,'TotalError'       );
    
end

function [d, x, y, A, I, Z] = preallocVectors (nCirc)

    %Initialize position vectors
    x = zeros(1, nCirc);
    y = zeros(1, nCirc);

    if nCirc==2
        d = 0;
        I = 0;    
        A = zeros(1,2);
        Z = zeros(1,3);

    else %nCirc==3
        d = zeros(1,3);
        I = zeros(1,4);
        A = zeros(1,3);
        Z = zeros(1,7);
    end
end

% *********************** END of Internal Function ********************

% *********************** Another Internal Function *******************
% *********************************************************************
function exitPos = ultraBars(inData, colors, rowNames, ...
    horizonBars , nextAxis)

% Input:
% inData: a matrix and vector of plotting data
% colors: colors for each unique value of inData
% rowNames: a cell array for name for each row
% horizonBar: specifies whether to add a horizontal bar to each plot or not
 
% get the number of unique numbers and remove the 0 which means not plot
uniqueVars = unique(inData) ;
uniqueVars(uniqueVars == 0) = [] ;

% create the legend variable if the inData is a row vector
if size(inData,1) == 1
    lgdVar = split( num2str(uniqueVars) )';
end

% get the number of colours to plot
if ~exist('colors','var') || isempty(colors)
    colors = rand(length(uniqueVars), 3);
end

% check the are is a row name for each row in inData
if size(inData,1) > 1 % only valid for matrix data
    if size(inData,1) ~= size(rowNames,1)
        error('row names should contain a name for each row in plotting data')
    end
end

% check if the orizontal bars are required on the plot
if ~exist('horizonBars','var')
    horizonBars = false;
end

% check for color errors 
if length(uniqueVars) > size(colors,1) 
    error('A color must be specified for each value in Data')
end

% plot the figure for row vectors 
% I have separated the two peaces of code because one does not need a box
% around the data
if size(inData,1) == 1
    % figure()
    % axes('position',[0.1300,0.11,0.7,0.04]);
    % set the bar width for large case
    if size(inData,2) < 500
        barwidth = 0.9 ;
    else
        barwidth = 1 ;
    end
    
    % now plot the data
    for ii = 1:length(uniqueVars)
        % get the data for that values and plot it in that color
        plotData = double(ismember(inData,uniqueVars(ii) )) ;
        bar(plotData,'FaceColor',colors(ii,:),...
            'EdgeColor',[1 1 1] ,'BarWidth',barwidth) ;
        hold on
    end
    set(gca,'GridColor',[1 1 1], 'XLim', [0.5 size(inData,2)+0.5 ], ...
        'XColor',[1 1 1] ,'YColor',[1 1 1],'XTick',[] , 'YTick',[],...
        'FontWeight','bold')
    % add a legend to the bar using the legendflex function
   
% if the is more thhan one row in the dataset
else
    % make a plot of multiple bar chart of top of each other
    % add the subtype plots to the heatmap using a loop
    % initialise some variables
    global plotTime
    figure(plotTime*100)
    % set(gcf,'position',[100,50,800,600])
        
    yInitial = 0.05; yPosSaved = yInitial;
    ySize = 0.02; increaseby = 0.022; % 0.44 0.4
    % change the size the bar if plotTime == 3
    if plotTime == 3
        ySize = 0.05; increaseby = 0.055;
    end
    xEndPos = 0.7 ; % 0.7750
    % loop over the rows and ascend by column
    for jj = 1:size(inData,1) % begin plotting from the bottem
        % define the next axis for the follow up plots
        if ~exist('nextAxis','var') || isempty(nextAxis)
            axes('position',[0.1300,yInitial,xEndPos,ySize]);
        elseif exist('nextAxis','var') && jj > 1
            axes('position',[0.1300,yInitial,xEndPos,ySize]);   
        else
            axes('position',nextAxis);
            yInitial = nextAxis(2) ;
            ySize = nextAxis(4) ; xEndPos = nextAxis(3) ;
            yPosSaved = yInitial ;
        end         
        for ii = 1:numel(uniqueVars) 
            plotData = double(ismember(inData(jj,:),uniqueVars(ii) )) ;
            bar(plotData,'FaceColor',colors(ii,:),'EdgeColor',[1 1 1] ,...
                'BarWidth',0.9) ;   
            hold on
            
            % add the name of the genes to the left of heatmap
            if exist('rowNames','var')
                dim = [0.02 yInitial 0.11 increaseby];
                annotation('textbox',dim,'String',rowNames{jj},...
                'FitBoxToText','on','FontSize',10,'EdgeColor','none',...
                    'HorizontalAlignment','right','FontWeight','bold',...
                    'VerticalAlignment','middle');
            end
        end
        % change the plot properties
        set(gca,'GridColor',[1 1 1], 'XLim', [0.5  size(inData,2)+0.5],...
            'XColor',[1 1 1] ,'YColor',[1 1 1],'YTickLabel',[],...
            'XTickLabel',[],'FontWeight','bold','YTick',[],'XTick',[])
        % increase the value to change the colors and plot positions
        yInitial = yInitial + increaseby;
    end
    % add a grey box to the plot usign the annotation function
    if plotTime ~= 3 % dont add the box if plotTime is == 3
    dim = [0.1300, yPosSaved, xEndPos, increaseby*size(inData,1)];
    annotation('rectangle',dim ,'Color',[0.5, 0.5, 0.5])
    end
    hold off 
   
    % plot the horizontal bars if they are required
    % prellocate the bar data size
    barhData = zeros(size(inData,1),numel(uniqueVars)) ;
    if horizonBars == true
        axes('position',[xEndPos+0.137, yPosSaved, 0.12, ...
            increaseby*size(inData,1)+0.001 ]);
        for kk = 1:numel(uniqueVars) 
            barhData(:,kk) = sum(inData == uniqueVars(kk),2) ;   
        end
        bar1 = barh(barhData,'stacked','BarWidth',0.85) ;
        % make sure there are no colors and spaces between the axis and the
        % first and last bar
        set(gca,'GridColor',[1 1 1], 'YLim', [0.5 size(inData,1)+0.52 ], ...
            'XColor',[0.5 0.5 0.5] ,'YColor',[1 1 1],'FontSize',10,...
            'YTick',[],'FontWeight','bold','Box','off','TickDir', 'out',...
            'LineWidth',1,'XAxisLocation','origin')
        
        % annoate the bar graph
        for ii = 1:size(barhData,2)
            set(bar1(ii),'FaceColor',colors(ii,:))
        end
        
    end
    exitPos = [0.1300,yInitial,xEndPos, ySize] ;
end

end

% *********************************************************************
% ********************** The first internal function ******************

function createLegendInternal(yPoint, xStart, legendLabels , plotColors,...
    myLgdTitle , fontSizes ,rectAndTextBox)

% specificy the y values starts and mode of actions for the drugs
% yPoint = 0.830 ; xStart = 0.1023 ;
xStartText = xStart + 0.01 ;
yPointTitle = yPoint + 0.03 ;

% specify the font size to be used in the plot
if ~exist('fontSizes','var')
    fontSizes = [10, 12] ;
end

% specifiy the rectangle and text length
if ~exist('rectAndTextBox','var')
    rectAndTextBox = [0.018 ,0.12] ;
end

% check for errors
if ~isnumeric(yPoint) || ~isnumeric(xStart)
    error('Both yPoint and xStarts should be numeric values')
elseif yPoint > 1 || xStart > 1
    error('Both yPoint and xStarts should be less than 1')
elseif ~isnumeric(plotColors)
    error('plot Color should be numeric')
end

if size(plotColors,1) ~= length(legendLabels)
    error('There should be a color for each legend names')
end

if iscategorical( legendLabels)
    legendLabels = categories(legendLabels);
end

for ii = 1:length(legendLabels)
    % add the legend color
    annotation('rectangle',[xStart yPoint rectAndTextBox(1) 0.023],...
        'EdgeColor', plotColors(ii,:), ...
        'FaceColor', plotColors(ii,:));
    
    % add the legend text
    annotation('textbox',[xStartText yPoint rectAndTextBox(2) 0.0230],...
        'String',legendLabels{ii},'FontSize',fontSizes(1),...
        'FontName','Helvetica Neue','FitBoxToText','off',...
        'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1] , ...
        'VerticalAlignment','middle','FontWeight','normal')
    
    % move the y point down
    yPoint = yPoint - 0.03 ;
end

% add the title
annotation('textbox',[xStart yPointTitle rectAndTextBox(2) 0.0230],...
    'String', myLgdTitle,'FontSize',fontSizes(2),...
    'FontName','Helvetica Neue','FitBoxToText','off',...
    'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1] , ...
    'VerticalAlignment','middle','FontWeight','bold',...
    'HorizontalAlignment','left');

end

% *********************** END of Internal Function ********************

% Transpose any table with this function 
%
% Input = any table
%
% Output = a transpose of the input table

function mytable = transposeTable(in_table)

myArray = table2cell(in_table(:,2:end) );
myArray = cell2table(myArray'); 
var_names = cellstr( table2cell(in_table(:,1)) );
var_names = matlab.lang.makeValidName(var_names) ;
var_names = var_names';

myArray.Properties.VariableNames = var_names ;

% S = {'my.Name','my_Name','my_Name'};
% validValues = matlab.lang.makeValidName(S)
% validUniqueValues = matlab.lang.makeUniqueStrings(validValues,{},...
%     namelengthmax)
  
row_names = in_table.Properties.VariableNames(2:end); 
row_names = cell2table(row_names');
mytable = [row_names , myArray ] ;
mytable.Properties.VariableNames(1,1) = ...
    in_table.Properties.VariableNames(1,1);

clear myArray var_names row_names ii str expression replace newStr 

end

% *********************** END of Internal Function ********************

% ****************************************************************

% ======================= another function =========================
function pSuperScript = convertPValue2SuperScript(p)
    % converts the number to scientific superscript for printing on a
    % figure
    pS = num2str(p) ;

    % get the first number
    firstNumbers = extractBefore(pS,'e') ;

    % check if there is a decimal place. then only get the first 4 numbers
    if contains( firstNumbers  ,'.')
        firstNumbers = firstNumbers(1:4) ;
    end
    
    % get the correctly formated p value
    pSuperScript = sprintf('%s x 10^{%d}', firstNumbers, ...
        str2double(extractAfter(pS, 'e') )) ;
    
    % if the p value is large
    if p > 0.0001
       pSuperScript = sprintf('%0.4f', p) ;
    elseif p == 0
         pSuperScript = sprintf('< 1 x 10^{%d}', -300) ;
    end

end

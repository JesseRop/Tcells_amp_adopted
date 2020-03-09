tabPanel(
  title = "About",
  value = "about",

  mainPanel(width = 12,
    
    HTML(
"
<h1>Welcome</h1>
<p>This site provides a platform for deeper interogation of sc-RNA outputs from the
<a href='https://www.gla.ac.uk/researchinstitutes/iii/'>Institute of Infection, Immunity & Inflammation (III), University of Glasgow</a>.
<p>


<h1>Respective publications</h1>
<p>
Please cite our paper if you use this data in your work:
</p>

<div class='citation'>
<h4>Beta2 integrins differentially regulate γδ T cell subset thymic development and peripheral maintenance</h4>
<p class='citation-authors'>C. L. McIntyre, L. Monin, T. D. Otto, C. S. Goodyear, A. C. Hayday, and V. L. Morrison.
</p>

<h4>Abstract</h4>
<p>γδ T cells reside predominantly at barrier sites and play essential roles in immune protection
against infection and cancer. Despite recent advances in the development of γδ T cell
immunotherapy, our understanding of the basic biology of these cells, including how their
numbers are regulated in vivo, remains poor. This is particularly true for tissue-resident γδ T cells.
We have identified the β 2 family of integrins as novel regulators of γδ T cells. β 2 integrin-deficient
mice displayed a striking increase in numbers of IL-17-producing Vγ6Vδ1 + γδ T cells in the lungs
and uterus, as well as circulation. Thymic development of this population was normal. However,
single cell RNA sequencing revealed the enrichment of genes associated with T cell survival and
proliferation specifically in β 2 integrin-deficient IL-17 + cells compared to their WT counterparts.
Indeed, β 2 integrin-deficient Vγ6 + cells from the lungs showed enhanced survival ex vivo,
suggesting that increased survival contributes to the accumulation of these cells in β 2 integrin-
deficient tissues. Furthermore, our data revealed an unexpected role for β 2 integrins in promoting
the thymic development of the IFNγ-producing CD27 + Vγ4 + γδ T cell subset. Together, our data
reveal that β 2 integrins play important roles in maintaining γδ T cell homeostasis, particularly in
regulating Vγ6 + cell numbers in mucosal tissues by controlling survival and Vγ4 + cell development
in the thymus. Our study indicates new and unprecedented mechanisms of control for γδ T cell
subsets.
</p>

<div class='text-center'>
<button id='button-data' type='button' class='btn btn-primary btn-lg'>Explore sc-RNA seq results for this work</button>
</div>
" 
# <h1>Contact us</h1>
# <p>
# Please <a href='tdo'>contact us</a>
# if you have any questions or comments on the analysis and results.
# </p>
# 
# <img src='synovial_pipeline.png' style='width:100%;'>
# </div>
# 
# <h1>Download the data</h1>
# <p>Download the complete data from ImmPort:</p>
# <ul>
# <li><a href='http://www.immport.org/immport-open/public/study/study/displayStudyDetail/SDY998'>SDY998:
# AMP Rheumatoid Arthritis Phase 1</a></li>
# </ul>
# 
# <h1>Get the code</h1>
# <p>Get the code for data analysis:
# </p>
# <ul>
# <li><a href='https://github.com/immunogenomics/amp_phase1_ra'>github.com/immunogenomics/amp_phase1_ra</a></li>
# </ul>
# <p>Get the code for this website:
# </p>
# <ul>
# <li><a href='https://github.com/immunogenomics/amp_phase1_ra_viewer'>github.com/immunogenomics/amp_phase1_ra_viewer</a></li>
# </ul>
# 
# 
# "
    )

    # h3("Accelerating Medicines Partnerships (AMP)"),
    # p(
    #   "The",
    #   a("Accelerating Medicines Partnership (AMP)",
    #     href = "https://www.nih.gov/research-training/accelerating-medicines-partnership-amp"),
    #   " is a public-private partnership between the National Institutes of",
    #   " Health (NIH), the U.S. Food and Drug Administration (FDA), 10",
    #   " biopharmaceutical companies and multiple non-profit organizations",
    #   " to transform the current model for developing new diagnostics and",
    #   " treatments by jointly identifying and validating promising",
    #   " biological targets for therapeutics. The ultimate goal is to",
    #   " increase the number of new diagnostics and therapies for patients",
    #   " and reduce the time and cost of developing them."
    # ),
    # 
    # h3("AMP RA Phase I"),
    # p(
    #   "Detecting distinct cellular subsets in tissues affected by rheumatoid arthritis is the key to deciphering pathogenesis in RA.",
    #   "We applied a multi-modal high dimensional strategy", 
    #   "including single-cell RNA-seq, mass cytometry, bulk RNA-seq, and flow cytometry",
    #   "to synovial tissue samples from 51 individuals with RA and osteoarthritis.",
    #   "Using an integrative strategy that uses canonical correlational analysis,", 
    #   "we are able to integrate across these data sets to define cellular populations that are robust.",
    #   "Evidence of these populations are seen across the different data modalities.",
    #   "This website supports the results of our AMP RA Phase I paper",
    #   a("(preprint version).",
    #     href = "https://www.biorxiv.org/content/early/2018/06/20/351130"),
    #     "Welcome to read the paper to know more details."
    # ),

    # h2("Disclaimer"),
    # p(
    #   "Data presented on this page is from Phase 1 of the AMP partership."
    #   # " Currently, this is private data meant to be shared internally,",
    #   # " only with consortium members."
    # ),
    # p(
    #    strong(
    #     "Sharing any data from this site with anyone outside of the",
    #    " AMP partnership is prohibited."
    #   )
    # ),
    # p(
    #   "This website is an experiment in providing early access to",
    #   " preliminary data analysis results. The content of this site is",
    #   " subject to change at any time without notice. We hope that you",
    #   " find it useful, but we provide it 'as is' without warranty of",
    #   " any kind, express or implied."
    # ),
    
    # h3("Contact"),
    # p(
    #   "This site is maintained by", 
    #   a("Kamil Slowikowski", href = "mailto:kslowikowski@fas.harvard.edu"),
    #   "and",
    #   a("Fan Zhang.", href = "mailto:fanzhang@broadinstitute.org"),
    #   "Please contact us if you have any questions, requests, or comments",
    #   " on the analysis and results."
    # ),

  ) # mainPanel
  
) # tabPanel

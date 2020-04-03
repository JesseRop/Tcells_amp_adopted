tabPanel("About",
         tags$head(tags$script(HTML('
        var fakeClick = function(tabName) {
          var dropdownList = document.getElementsByTagName("a");
          for (var i = 0; i < dropdownList.length; i++) {
            var link = dropdownList[i];
            if(link.getAttribute("data-value") == tabName) {
              link.click();
            };
          }
        };
      '))),
         mainPanel(width = 12,
                   HTML("
<h1>Welcome</h1>
<p>
Here, you can explore the scRNAseq dataset from wild type (WT) and β2 integrin-deficient (CD18 Knock-out[KO]) γδ T cells isolated from mice lungs.
</p>
<p>
Further details can be found in a our paper summarized below which you should cite as relevant
</p>

<div class='citation'>
<h4><a href='https://www.google.com' target='_blank'>Beta2 integrins differentially regulate γδ T cell subset thymic development and peripheral maintenance</a></h4>
<p class='citation-authors'>C. L. McIntyre, L. Monin, T. D. Otto, C. S. Goodyear, A. C. Hayday, and V. L. Morrison.
</p>
</div>


<div>

<h4>Paper significance statement</h4>
<p>γδ T cells reside in barrier tissues and provide immune protection against infection and cancer. Their anti-tumor potential has
led to recent advances in the development of γδ T cell immunotherapy. However, our understanding of the basic biology of these 
cells, including what molecules and pathways control their maintenance within barrier tissues, remains poor. We demonstrate that
β2 integrin adhesion molecules play a major role in regulating γδ T cell subset numbers during homeostasis: the loss of β2 integrin
expression results in a striking increase in IL17-producing γδ T cells in the lungs and uterus due to enhanced survival. These 
findings illustrate a novel mechanism of γδ T cell regulation that may have significant implications for immunotherapy development.
</p>
</div>

<div class='citation'>
<p> The results are displayed in 2 sections
<ul>
<li><a onclick = 'fakeClick(\"Differential expression\")'> Differential expression</a> - This allow the visualization of gene expression in CD18 KO cells to that of Wildtype cells using
parameter specified in the paper (the first 15 Principal components and 0.15 as the clustering resolution). </li>
<li><a onclick = 'fakeClick(\"Cluster adjustment\")'> Clustering</a> - This allows advanced analysis by increasing the resolution to further refine the clustering of cells. Increased resolution results in
in more clusters which may be sub-population of previous clusters.</li>
</ul>
</p> 
</div>

<div>
<h4>Key results</h4>
<img src='tcell_cd18_intro_image.png' alt='Results' style='width:100%;'>
</div>

"))
         
)


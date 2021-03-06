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
Here, you can explore the single-cell RNA sequencing (scRNA-Seq) results from wild type (WT) and β2 integrin-deficient (CD18 Knock-out[KO]) gamma delta (γδ) T cells isolated from mice lungs.
</p>
<p>
The analysis was carried out using the Seurat pipeline for differential expression. Further details can be found in our paper summarized below which you should cite if you use our data.
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
<li><a onclick = 'fakeClick(\"Differential expression\")'> Differential expression (DE)</a> - Visualization of gene expression in CD18-KO cells compared to that in WT cells using
parameters specified in our <a href='https://www.google.com' target='_blank'> paper </a> (i.e the first 15 Principal components and 0.15 as the clustering resolution). </li>
<li><a onclick = 'fakeClick(\"Cluster adjustment\")'> Clustering</a> - Allows flexibility in increasing the number of clusters by increasing the resolution. 
This results in the sub-division of clusters characterized in our <a href='https://www.google.com' target='_blank'> paper </a> into sub-populations.</li>
</ul>
</p> 
</div>

<div style='border:thin solid black ; padding:30px ; margin:30px'>
<figure>  
<h4>Key results</h4>
<img src='tcell_cd18_intro_image_clearer_labelled.png' alt='Results' style='width:100%;'>
<figcaption><b>A:</b> UMAP of the labelled  populations of γδ T cells obtained using 15 principal components and 0.15 as the clustering resolution in Seurat. <b>B:</b> The proportion 
of the different populations of γδ T cells in the wildtype replicates (WT1 and WT2) and the CD18 Knockout (KO1).</figcaption>
</figure>
</div>

<footer align = 'center' style = 'position:absolute; bottom:80; width:95%; height:20px; padding: 10px; z-index: 200; background-color:PowderBlue;'> <hr></hr> <p>This web-app was developed by <a href='https://https://www.gla.ac.uk/researchinstitutes/iii/staff/?webapp=staffcontact&action=person&id=4edfe8ec8a95' target='_blank'> Jesse Rop </a> as a <a href='https://https://www.ukri.org/research/global-challenges-research-fund/' target='_blank'> GCRF </a> exchange program project under the supervision of <a href='https://https://www.gla.ac.uk/researchinstitutes/iii/staff/thomasdanotto/' target='_blank'> Dr Thomas Otto </a> and with the assistance of <a href='https://https://www.gla.ac.uk/researchinstitutes/iii/staff/vickymorrison/' target='_blank'> Dr Vicky Morrison </a> and <a href='https://https://www.gla.ac.uk/researchinstitutes/iii/staff/?webapp=staffcontact&action=person&id=49d8ebe38396' target='_blank'> Scott Arkinson </a>. Further acknowledgements on the generation of the data hosted in the app can be found in our <a href='https://www.google.com' target='_blank'> paper </a>.</p> <p> The code used to develop the app can be found <a href='https://github.com/JesseRop/Tcells_amp_adopted/blob/master/introduction_page.R' target='_blank'> here </a>. For any queries or suggestions kindly <a href='mailto:ThomasDan.Otto@glasgow.ac.uk' target='_blank'> contact us</a>.</p>
</footer>"))
         
)

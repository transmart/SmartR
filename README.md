[![Build Status](https://travis-ci.org/transmart/SmartR.svg?branch=master)](https://travis-ci.org/transmart/SmartR)

###SmartR is a grails plugin seeking to improve the visual analytics of the [tranSMART platform](https://github.com/transmart/transmartApp) by using recent web technologies such as [d3](http://d3js.org/).

[![Heatmap Example](https://i.imgur.com/WGFV2kD.png)](https://youtu.be/Gg0AdYt77Cs)

#### Installation
##### tranSMART 16.2, eTRIKS v3 and younger builds
SmartR will be included by default. No modification necessary.

##### tranSMART 1.2.5 and 16.1
[Just modify your build configuration](https://github.com/transmart/transmartApp/commit/a781506f06d6ffb38cba307ed13879cb35e22056).

##### tranSMART 1.2.4 and older
[Modify your build configuration](https://github.com/transmart/transmartApp/commit/a781506f06d6ffb38cba307ed13879cb35e22056).
And add the following code to transmartApp/web-app/js/datasetExplorer/datasetExplorer.js near [this location:](https://github.com/transmart/transmartApp/blob/release-1.2.4/web-app/js/datasetExplorer/datasetExplorer.js#L782)
```
loadPlugin('smartR', "/SmartR/loadScripts", function () {
    resultsTabPanel.insert(4, smartRPanel);
})
```

#### Requirements
- SmartR requires the following R packages:
  - [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)
  - [gplots](https://cran.r-project.org/web/packages/gplots/index.html)
  - [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
  - [WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html)
  - [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
  
#### Releases
- https://github.com/transmart/SmartR/releases

#### Wiki
- https://wiki.transmartfoundation.org/display/transmartwiki/SmartR

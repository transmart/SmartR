###SmartR is a grails plugin seeking to improve the visual analytics of the [tranSMART platform](https://github.com/transmart/transmartApp) by using recent web technologies such as [d3](http://d3js.org/).

[![Heatmap Example](https://i.imgur.com/WGFV2kD.png)](https://youtu.be/Gg0AdYt77Cs)

#### Installation
##### tranSMART 16.2, eTRIKS v3 and younger builds
SmartR will be included by default. No modification necessary.

##### custom builds
1. Add ```runtime ':smart-r:1.2-STABLE-SNAPSHOT'``` to **BuildConfig.groovy** in the **transmartApp** source code.
2. Compile a WAR file via ```grails war``` for deployment.

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

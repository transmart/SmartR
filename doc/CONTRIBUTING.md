## SETUP (The following steps HIGHLY depend on your environment. Feel free to skip those if you know what you are doing or have a better tutorial)

0. Preparation
    - Download [Grails 2.3.11](https://grails.org/download.html)
    - Download [Oracle JDK7](http://www.oracle.com/technetwork/java/javase/downloads/jdk7-downloads-1880260.html) or OpenJDK 7
        - Version 8+ is not supported at the moment.
    - Download the [tranSMART configuration template](https://github.com/transmart/transmart-data/blob/master/config/Config-template.groovy)
    - Download the [tranSMART data source configuration](https://github.com/transmart/transmart-data/blob/master/config/DataSource.groovy.php)
    - Clone the [tranSMART repository](https://github.com/transmart/transmartApp)
    - Clone the [SmartR repository](https://github.com/transmart/SmartR)
    - Install the [necessary R packages](https://github.com/transmart/SmartR#requirements)

1. Add JAVA_HOME, GRAILS_HOME, and the grails binaries to your environment
In an unix-like environment you could add something like this to your .bashrc
    ```
    export JAVA_HOME=/Library/Java/JavaVirtualMachines/jdk1.7.0_80.jdk/Contents/Home/
    export GRAILS_HOME=~/.grails/grails-2.3.11/
    export PATH=$GRAILS_HOME/bin:$PATH
    ```

2. If you use Oracle JDK7 you must update your cacerts file, with a newer version (i.e. from Oracle JDK8)
    ```
    sudo cp /Library/Java/JavaVirtualMachines/jdk1.8.0_112.jdk/Contents/Home/jre/lib/security/cacerts /Library/Java/JavaVirtualMachines/jdk1.7.0_80.jdk/Contents/Home/jre/lib/security/cacerts
    ```

3. Prepare your .grails folder
    ```
    mkdir ~/.grails
    mkdir ~/.grails/transmartConfig
    mv Config-template.groovy ~/.grails/transmartConfig/Config.groovy
    mv DataSource.groovy.php ~/.grails/transmartConfig/DataSource.groovy
    ```

4. Link tranSMART to your SmartR directory
    - Search this line `runtime ':resources:1.2.1'` and change it to `runtime ':resources:1.2.14'`
    - Remove this line `runtime ':smart-r:16.2-SNAPSHOT'`
    - Add this line to the end of the file `grails.plugin.location.SmartR="../relative/path/to/SmartR"`

5. Prepare your Configuration files
    - Add this line to the end of Config.groovy `org.transmart.configFine = true`
    - Modify parameters if necessary
    - Remove the php code from DataSource.groovy and edit it as needed

6. (optional) If you want to use a remote DB, you must establish an connection to it now. This step depends on your settings in DataSource.groovy
    ```
    ssh -L 3080:localhost:5432 foobar@10.11.12.13 -p 8022
    ```

7. Start Rserve
    ```
    R
    require(RServe)
    Rserve(args="--no-save")
    ```

8. Run tranSMART
    ```
    cd /wherever/transmartApp
    grails run-app
    ```
    NOTE: This will crash the first time, probably because of some circular dependency. Run it again and it should work.

## SmartR development

1. Coding Guidelines
    - The code must pass all checks of jshint, with the ruleset defined in .jshintrc

2. Tests
    - Your code must pass all tests in order to be accepted

## SmartR Components & Workflow Structure (Example Boxplot)

Note: Some parts of a workflow are identified by pattern matching. When you define a new workflow you should think of a name like "MyCoolWorkflow" and apply that name as shown in the Boxplot example.

#### Register your workflow in `grails-app/conf/SmartRResources.groovy`
```
resource url: [plugin: 'smart-r', dir: 'js/smartR/_angular/controllers', file: 'boxplot.js']
resource url: [plugin: 'smart-r', dir: 'js/smartR/_angular/viz', file: 'plotlyBoxplot.js']
```

#### Code your vizualisation in `web-app/js/smartR/_angular/viz/plotlyBoxplot.js`

This is a standard angular directive that watches the model where the R results will be stored. The vizualization code is stored within the link function.

#### Add your workflow controller in `web-app/js/smartR/_angular/controllers/boxplot.js`

This is the controller that "manages" the workflow elements. It is recommended to stay close to the already existing workflows for consistent behaviour.

#### Apply your css in `web-app/css/boxplot.css`

Just dome CSS for the code in plotlyBoxplot.js

#### Define the workflow layout in `grails-app/views/layouts/_boxplot.gsp`

There are several pre-defined components available to define your workflow layout.
Here are some of the most important ones. While there are some other useful ones this will cover your most basic needs.
For more details inspect the directive code directly.

```
<tab-container>
Wraps several <workflow-tab> and creates the "tabbed view".
```

```
<workflow-tab>
A single tab in a "tabbed view".

Attributes:
tab-name: Name of the tab
disabled: If true the tab greys out and is not selectable
```

```
<concept-box>
The place where you drag your nodes into.

Attributes:
concept-group: Must point to an object with two attributes, 'concepts'(Array) and 'valid'(Boolean), that is part of another object. See any other workflow for examples.
type: Can be "HD", "LD-numerical" or "LD-categorical"
min: minimum amount of nodes for this box
max: maximum amount of nodes for this box (-1 if no limit)
label: Label of the box
tooltip: Tooltip of the box
```

```
<fetch-button>
This button loads all data specified in the <concept-box>s and loads them into the RSession.

Attributes:
concept-map: Points to the parent object of the object specified in concept-group in <concept-box>. See any other workflow for examples.
loaded: True if data have been fetched
running: True if fetching in progress
biomarkers: Points at the same object than <biomarker-selection>. This is for selecting only certain rows from HDD.
disabled: True if button should be greyed out
message: Message to show next to the button (errors, infos, ...)
allowed-cohorts: Can be [1], [2], or [1,2].
```

```
<run-button>
This button runs the run.R main() method and stores the result in the specified model.

Attributes:
running: True if R code is running
store-results-in: The model where the script results will go
button-name: Name of the button
filename: If you are generating a file instead of returning json directly from the R script, then the name of that file should be specified here. (e.g. heatmap.json)
params: Arguments that you would like to call the main() method with
```

#### Add the R code in `web-app/HeimScripts/boxplot/run.R`

The R code here will be executed by clicking the `<run-button>`.
Every .R file that you execute from SmartR must have a main() function, which is the entry point for each call to that script.
In all .R files executed by SmartR there will be two global variables present: `loaded_variables` and `fetch_params`.

`loaded_variables` contains the fetched data and is a named list of data frames, where the names follow this conventions:
{CONCEPTS}_n{NODEID}_s{SUBSET}, where CONCEPTS matches the name given in the controller, NODEID the index of the node within the concept box, and SUBSET the subset number.
An example would look like this: hddConcepts_n4_s2

`fetch_params` contains several meta data for every node.

The main method must either return valid JSON (e.g. `toJSON(someList)`) or a file if your workflow makes use of high dimensional data (e.g. `write(toJSON(someList), file="boxplot.json")`)

Tip: `save(loaded_variables, file="~/loaded_variables.Rda")` and `load("loaded_variables.Rda")` are really handy if you want to develope the R code seperately.

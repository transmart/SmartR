<r:require modules="smartR_heatmap"/>
<r:layoutResources disposition="defer"/>

<div style="width: 98%">

    <p>
        Start analysis with <b>Fetch</b>, <b>Preprocess</b>  and then <b>Run Analysis</b>. Make sure that you run Summary
    Statistic prior running the Heatmap.
    </p>

    <div id="heim-tabs" style="margin-top: 25px;">

        <ul>
            <li><a href="#fragment-load"><span>Fetch</span></a></li>
            <li><a href="#fragment-preprocess"><span>Preprocess</span></a></li>
            <li><a href="#fragment-run"><span>Run</span></a></li>
        </ul>

        %{--========================================================================================================--}%
        %{--Load Data--}%
        %{--========================================================================================================--}%
        <div id="fragment-load">

            <div class="sr-input-area">

                %{--High dimension dropzone--}%
                <div class="heim-input-field heim-dropzone sr-hd-input">
                    %{--High dimensional input--}%
                    <label>High dimensional <i class="ui-icon ui-icon-info" title="Select high dimensional data node(s) from the Data
                        Set Explorer Tree and drag it into the box. The nodes needs to be from the same platform.">
                    </i></label>

                    <br><br>

                    <div id='divIndependentVariable'
                         style="border:1px solid #666; height: 100px; padding: 10px;"></div>

                    <div style="margin-top: 10px; text-align: right;">
                        <button type="button" id="heim-btn-clear">Clear</button>
                    </div>
                </div>

                %{--TODO: Support numerical nodes--}%
                %{--Low dimension dropzone (numerical)--}%
                <div class="heim-input-field heim-dropzone sr-low-input">
                    %{--Numerical node input--}%
                    <label>Numerical (optional) <i class="ui-icon ui-icon-info" title="Select numerical data node from the Data
                            Set Explorer Tree and drag it into the box.">
                    </i></label>
                    <br><br>

                    <div id='divNumVariable'
                         style="border:1px solid #CCC; height: 100px; padding: 10px;"></div>

                    <div style="margin-top: 10px; text-align: right;">
                        <button type="button" id="heim-btn-clear-2" disabled>Clear</button>
                    </div>
                </div>

                %{--TODO: Support categorical  nodes--}%
                %{--Low dimension dropzone (categorical)--}%
                <div class="heim-input-field heim-dropzone sr-low-input">
                    %{--High dimensional input--}%
                    <label>Categorical (optional) <i class="ui-icon ui-icon-info" title="Select categorical data node from the Data
                                Set Explorer Tree and drag it into the box.">
                    </i></label>
                    <br><br>

                    <div id='divCategoricalVariable'
                         style="border:1px solid #CCC; height: 100px; padding: 10px;"></div>

                    <div style="margin-top: 10px; text-align: right;">
                        <button type="button" id="heim-btn-clear-3" disabled>Clear</button>
                    </div>
                </div>
            </div>

            %{--TODO: Support gene filtering--}%
            <div class="sr-fetch-params-area">
                %{--Select identifier--}%
                <div class="heim-input-field heim-autocomplete" style="clear: both">
                    <label>Select a Gene/Pathway/mirID/UniProtID:</label>
                    <input id="heim-input-txt-identifiers" disabled>
                </div>
            </div>

            %{--tool buttons--}%
            <div style="margin-top: 10px; text-align: right;">
                <button id="heim-btn-fetch-data">Fetch Data</button>
            </div>

            <div id="heim-fetch-data-output"></div>

        </div>

        %{--========================================================================================================--}%
        %{--Preprocess--}%
        %{--========================================================================================================--}%
        <div id="fragment-preprocess">

            %{--Aggregate Probes--}%
            <div class="heim-input-field">
                <input type="checkbox" id="chkAggregateProbes" name="chkAggregateProbes">
                <span>Aggregate probes</span>
            </div>

            %{--tool buttons--}%
            <div style="margin-top: 10px; text-align: right;">
                <button id="heim-btn-preprocess-heatmap">Preprocess</button>
            </div>

            %{--Preprocess Output--}%
            <div id="heim-preprocess-output"></div>
        </div>

        %{--========================================================================================================--}%
        %{--Run Analysis--}%
        %{--========================================================================================================--}%
        <div id="fragment-run">

            %{--Number of max row to display--}%
            <div class="heim-input-field heim-input-number">
                <label>Number of max row to display</label>
                <input type="text" id="txtMaxRow" value="100">
            </div>


            %{--TODO :Not relevant for now. Might be in future sprint--}%

            %{--Group Subject--}%
            %{--<div class="heim-input-field">--}%
            %{--<input type="checkbox" id="chkGroupSubject" name="chkGroupSubject">--}%
            %{--<span>Group by subject instead of node (only applicable when multiple nodes are selected).</span>--}%
            %{--</div>--}%

            %{--Apply statistical methods--}%
            %{--<div class="heim-input-field heim-input-number">--}%
            %{--<label>Apply statistical methods</label>--}%
            %{--<select id="methodSelect">--}%
            %{--<option value="none">None</option>--}%
            %{--<option value="hierarchical-clustering">Hierarchical Clustering</option>--}%
            %{--<option value="k-means-clustering">K-Means Clustering</option>--}%
            %{--<option value="marker-selection">Marker Selection</option>--}%
            %{--</select>--}%
            %{--</div>--}%

            <hr style="margin-top: 20px;">

            %{--Options--}%
            <div id="clusteringOptionsDiv" hidden>
                %{--Number of clusters--}%
                <div class="heim-input-field heim-input-number" id="noOfClustersDiv">
                    <label>Number of clusters</label>
                    <input type="text" id="txtNoOfClusters" name="txtNoOfClusters" value="2">
                </div>

                %{--Number of markers--}%
                <div class="heim-input-field heim-input-number" id="noOfMarkersDiv">
                    <label>Number of markers</label>
                    <input type="text" id="txtNoOfMarkers" name="txtNoOfMarkers">
                </div>

                %{--apply row clustering--}%
                <div class="heim-input-field">
                    <input type="checkbox" id="chkApplyRowClustering" name="chkApplyRowClustering">
                    <span>Apply clustering for row.</span>
                </div>

                %{--apply column clustering --}%
                <div class="heim-input-field">
                    <input type="checkbox" id="chkApplyColumnClustering" name="chkApplyColumnClustering">
                    <span>Apply clustering for column.</span>
                </div>

            </div>

            %{--tool buttons--}%
            <div style="margin-top: 10px; text-align: right;">
                <button id="heim-btn-run-heatmap">Get Heatmap</button>
            </div>

            %{--result--}%

            <div id="heim-run-output"></div>

            %{--d3 js heatmap placeholder--}%
            <div id='visualization' class='text'>
                <div id="heatmap" class='text'></div>
            </div>
        </div>

    </div>
</div>

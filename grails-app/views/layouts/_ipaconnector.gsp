<script type="text/ng-template" id="ipaconnector">
<div ng-controller="IpaconnectorController as ipaconn">

    <tab-container>
        %{--========================================================================================================--}%
        %{-- Fetch Data --}%
        %{--========================================================================================================--}%
        <workflow-tab tab-name="Fetch Data" disabled="fetch.disabled">
            <concept-box
                concept-group="fetch.conceptBoxes.highDimensional"
                type="HD"
                min="1"
                max="-1"
                label="High Dimensional"
                tooltip="Select high dimensional data node(s) from the Data Set Explorer Tree and drag it into the box.
                The nodes needs to be from the same platform."
            >
            </concept-box>

            <hr class="sr-divider">
            <fetch-button
                    loaded="fetch.loaded"
                    running="fetch.running"
                    concept-map="fetch.conceptBoxes"
                    biomarkers="fetch.selectedBiomarkers"
                    show-summary-stats="false"
                    summary-data="fetch.scriptResults"
                    all-samples="common.totalSamples"
                    allowed-cohorts="[1,2]"
                    number-of-rows="common.numberOfRows"
            >
            </fetch-button>
        </workflow-tab>

        %{--========================================================================================================--}%
        %{--Run Analysis--}%
        %{--========================================================================================================--}%
        <workflow-tab tab-name="Run Analysis" disabled="runAnalysis.disabled">
            <div>
                NOTE: This version of the IPA Connector is a pilot and
                we would like to collect any feedback for future
                development. The IPA Connector requires access to a
                valid IPA license (Ingenuity Pathway Analysis or IPA
                is available from QIAGEN at qiagen.com) with the
                Export API enabled. If you have a valid license and do
                not have the Export API enabled please contact your
                QIAGEN Bioinformatics sales representative or QIAGEN
                Customer Support at:<br>
                <a href="mailto:AdvancedGenomicsSupport@qiagen.com?Subject=tranSMART%20IPA%20Connector">AdvancedGenomicsSupport@qiagen.com</a><br>
                U.S.: +1 866 464 3684<br>
                Germany: +49 (0) 341 3397 5301<br>
                Denmark: +45 80 82 0167
            </div>

            <form name="ipaconn.degParamForm">
                <div class="heim-input-field sr-input-area">
                    <div class="heim-input-field-sub" id="sr-multi-subset" >
                        <fieldset class="heim-radiogroup" id="sr-differential-exp-group" ng-disabled="subsets < 2">
                            <h3>Significance measure</h3>
                            <div><label>
                                <input type="radio" ng-model="runAnalysis.params.significanceMeasure" value="pval"> p-value
                            </label></div>
                            <div><label>
                                <input type="radio" ng-model="runAnalysis.params.significanceMeasure" value="adjpval"> Benjamini & Hochberg
                            </label></div>
                        </fieldset>
                    </div>
                </div>

                <div class="ipa-inputtext">
                    <label for="significanceCutoff">Significance cutoff</label>
                    <input type="number" min="0" max="1" step="any" id="significanceCutoff"
                           name="significanceCutoff" ng-model="runAnalysis.params.significanceCutoff"><br/>
                    <label for="fcCutoff">Fold-change cutoff</label>
                    <input type="number" min="0" step="any" id="fcCutoff"
                           name="fcCutoff" ng-model="runAnalysis.params.fcCutoff">
                </div>
            </form>

            <form name="ipaCredentials" id="ipa-credentials" class="ipa-inputtext">
                <label for="ipa-username">IPA username</label>
                <input type="text" id="ipa-username" ng-model="runAnalysis.ipaCredentials.username"/><br/>
                <label for="ipa-password">IPA password</label>
                <input type="password" id="ipa-password" ng-model="runAnalysis.ipaCredentials.password"/>
            </form>
            <div class="ipaSecurityWarning">
                Note: The IPA username and password will be sent as
                plain text to the tranSMART server.
                <span ng-show="! runAnalysis.ipaConnectionIsSecure"><span class="ipa-warning">WARNING</span>
                    Your connection does not seem to be
                    https-encrypted, this means an eavesdropper can
                    potentially read your IPA credentials. Only click
                    the <span class="ipa-code">Launch IPA</span> button if you
                    trust the entire route from your computer to the
                    tranSMART server (eg inside a company network).
                </span>
            </div>

            <hr class="sr-divider" />

            <run-button button-name="Calculate differential expression"
                        store-results-in="runAnalysis.scriptResults"
                        script-to-run="run"
                        arguments-to-use="runAnalysis.params"
                        filename="ipaconnector.json"
                        running="runAnalysis.running"
                        wait-message="Calculating differentially expressed genes, please wait"
            >
            </run-button>

            <ipa-api
                differentially-expressed="differentiallyExpressed"
                deg-params="runAnalysis.params"
                credentials="runAnalysis.ipaCredentials"
            >
            </ipa-api>
        </workflow-tab>

    </tab-container>
</div>
</script>

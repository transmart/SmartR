<!DOCTYPE html>
<html>
  <head>
    <title><g:if env="development">Grails Runtime Exception</g:if><g:else>Error</g:else></title>
    <meta name="layout" content="main">
    <g:if env="development"><link rel="stylesheet" href="${resource(dir: 'css', file: 'errors.css')}" type="text/css"></g:if>
  </head>
  <body>
    <g:logMsg>SmartR An error has occurred</g:logMsg>
    <g:logMsg>exception '${exception}'</g:logMsg>
    <g:if env="development">
      <g:renderException exception="${exception}" />
    </g:if>
    <g:else>
      <ul class="errors">
        <li>SmartR</li>
	<li>An error has occurred</li>
	<li><g:renderException exception="${exception}"/></li>
      </ul>
    </g:else>
  </body>
</html>

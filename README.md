# project_template

A skeleton template for analysis projects. To make full use of the reporting capabilities make sure Sphinx is installed locally. A workflow file and Buildcode have been added to publish a GitHub Pages website from the results of a project.

- Use this repo as a template  "Use this template" button. 
- Clone the project repo to your local computer.
- To activate the GitHub pages Sphinx static site, go to the Settings of the online project repo and activate GitHub pages from the Branch "gh-pages": every push to "main" will redeploy the Pages website
- Make html report locally by running "make html" in the docs/ folder of the local project repo for quick updates and checks. 
- Figures and Tables that are outputted from the scripts in the Scripts folder should be saved in the results folder; first intermediate and then in the final results section. If you want to publish them in the html report, save them in the docs/_static folder. References from the .rst files will refer to this folder.

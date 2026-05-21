This folder will describe the different measures we took to utilize the OligoAI model and data in our project.

The setup of the environment of OligoAI model is described in the [RunningModel.md](RunningModel.md) file.

We ran OligoAI on its own data, to get their predictions.

We also ran OligoAI on data from ASOptimizer study which they did not use.
To see how we adapted the data to run on OligoAI, please view the [ToOligoAIFormat.ipynb](../notebooks/data/ASOPtimizer/ToOligoAIFormat.ipynb), and how we ran the model here [OligoAI_with_ASOptimizer_data.ipynb](../notebooks/competitors/OligoAI/OligoAI_with_ASOptimizer_data.ipynb). To be able to run it we had to add context with our logic, utilizing the GRCh38 genome and a corresponding annotation file, available to view in the [cli.py](../src/tauso/cli.py) file. Notice there was some overlap between the ASOptimizer data and OligoAI data, 

We can view the unique data in the CompareData.ipynb notebook.


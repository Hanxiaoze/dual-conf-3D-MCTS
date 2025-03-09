# dual-conf-3D-MCTS
dual conformation 3D Monte Carlo Tree Search drug generation method with mini-GPT


<img width="1248" alt="TOC" src="https://github.com/user-attachments/assets/34dd3e54-7d01-49d5-9c6b-a3cd28a4361e" />




# Dependencies
------------
# Use envi file:
```
conda env create -f ./environment.yml -n dual-conf-3D-MCTS_env
```

# Main Environment dependencies:
```
# software
gnina         ## https://github.com/gnina/gnina
ADFR          ## https://ccsb.scripps.edu/adfr/
MPerformer    ## https://github.com/FanmengWang/MPerformer

# python
python >= 3.7                ## 
openbabel >= 3.1.1           ## conda install openbabel -c openbabel
rdkit >= 2022.03.2           ## conda install rdkit -c rdkit
func-timeout >= 4.3.5        ## pip install func-timeout


```
The GNINA we used was binary version built by the authors on Mar 6, 2021. It can be downloaded [here](https://drive.google.com/file/d/1m6Uf3ALlEnvgztEzZrcy7gVO4Ag7VI6-/view?usp=drive_link).

Create `MPerformer_test` conda env:
```bash
   conda create -n MPerformer_test python=3.9 pytorch=2.0 torchvision cudatoolkit=11.7 ase openbabel -c pytorch -c openbabel -c defaults -c conda-forge
   (download_file https://github.com/dptech-corp/Uni-Core/releases/download/0.0.3/unicore-0.0.1+cu117torch2.0.0-cp39-cp39-linux_x86_64.whl)
   conda run -n MPerformer_test pip install ./unicore-0.0.1+cu117torch2.0.0-cp39-cp39-linux_x86_64.whl
   conda run -n MPerformer_test pip install rdkit-pypi==2021.9.4
   conda run -n MPerformer_test pip install dpdata
   conda run -n MPerformer_test pip install torch==2.0 torchvision torchaudio
   conda run -n MPerformer_test pip install pandas
   conda run -n MPerformer_test pip install scikit-learn
   conda run -n MPerformer_test pip install numpy
```


Quick Start
------------
The script `dual_target_drug_prepare_traning_generating_analysis_workflow.py` is the main interface, it will automatically solve and creat Conda env and call other script to run and the whole workflow. The absolute paths of `MGLToolsPckgs` and `Vinna_bin_path` are need to change with your installing path. 

If you have solved the Conda env manually, you can just use the script `dual_target_drug_prepare_traning_generating_analysis_workflow_manu_solve_conda_env.py` to run the whole workflow.

Please download the [MPerformer-checkpoint](https://drive.google.com/file/d/1sHWm1xOy0I8_R50dPANfMUXoRQkoPCBJ/view?usp=drive_link) and place it to the fold `./MPerformer-master/weight`

In short, the workflow can be run by just two ways:

```bash
   python dual_target_drug_prepare_traning_generating_analysis_workflow.py
   python dual_target_drug_prepare_traning_generating_analysis_workflow_manu_solve_conda_env.py
   ```

The PDB files `target_protein_1.pdb` and `target_protein_2.pdb` are the 3D structures of the two target-proteins which are used as docking receptors, they should be changed to your target-proteins PDB files according to your researching project. The files `config_1.txt` and `config_2.txt` are two docking configuration files for target-protein-1 and target-protein-2 respectively, they definite the docking pocket sites of two target-proteins. 

The files `times.ttf`, `timesi.ttf`, `timesbi.ttf`, `timesbd.ttf` are the font files to set the matplotlib figure plotting font, the need to be placed to the right place of the corresponding Conda env (`~/anaconda3/envs/3D_Scaffold_test/lib/python3.8/site-packages/matplotlib/mpl-data/fonts/ttf/`).


Data for research reproducing
--------
1. [zip files for research reproducing of our paper of training and generating process](https://drive.google.com/drive/folders/1WnKOa9ul7w6HpIXFdSPEOsk43Ec5tr0i?usp=sharing)
2. [all files for research reproducing of our paper of MD simulations and MM-GBSA calculations](https://pan.baidu.com/s/18DbmaKho242RSxJ-cbmrPg?pwd=8888)

Reference
--------
1. Reference the single-target 3D-MCTS lib from: https: https://github.com/Brian-hongyan/3D-MCTS
2. Reference the gnina from: https://github.com/gnina/gnina
3. Reference the MPerformer lib from: https://github.com/FanmengWang/MPerformer
4. Reference the Uni-Core from: https://github.com/dptech-corp/Uni-Core/releases


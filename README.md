# dual-conf-3D-MCTS
dual conformation 3D Monte Carlo Tree Search drug generation method with mini-GPT


<img width="1248" alt="TOC" src="https://github.com/user-attachments/assets/34dd3e54-7d01-49d5-9c6b-a3cd28a4361e" />




# Dependencies
------------
# Use our provided envi file for easy create:
```
conda env create -f ./environment.yml -n dual-conf-3D-MCTS_env
conda env create -f ./mperformer_environment.yml -n mperformer_env
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
Please download the [MPerformer-checkpoint](https://drive.google.com/file/d/1sHWm1xOy0I8_R50dPANfMUXoRQkoPCBJ/view?usp=drive_link) and place it to the fold `./MPerformer-master/weight`



# Quick Start
------------
Train the mini-GPT/QWEN model for dual-conf-3D-MCTS use our provided file:
```
mini_QWEN_frag_recomd.ipynb
```
Then, you can run ```3D-MCTS-dual-conf_QWEN.py```

You need to in 3rd level sub-directory to run the script:
```
cd dual-conf-3D-MCTS/
conda activate dual-conf-3D-MCTS_env
mkdir -p first/second/run_dua_conf
cd ./first/second/run_dua_conf
cp ligand_A.sdf ligand_B.sdf conf_A.pdb conf_B.pdb conf_A_pocket.pdb conf_B_pocket.pdb ./   # copy your mol files to here
python ../../../3D-MCTS-dual-conf_QWEN.py --num_sims 100000 --ligand ./ligand_A.sdf --ligand_dual ./ligand_B.sdf --protein ./conf_A.pdb --protein_dual ./conf_B.pdb --pocket ./conf_A_pocket.pdb --pocket_dual ./conf_B_pocket.pdb --score -7 --score_dual 0 --start 1 --qed 0.3 --processor 30
```

An example, for the paper data reproduce, dual conformation:
```
cd dual-conf-3D-MCTS/SARS-CoV-2/Mpro_A_B_conformation_generate_analysis/run_dual_QWEN_test/
python ../../../3D-MCTS-dual-conf_QWEN.py --num_sims 100000 --ligand ../ligand_A_s.sdf --ligand_dual ../ligand_B_s.sdf --protein ../conf_A_919.pdb --protein_dual ../conf_B_4240.pdb --pocket ../conf_A_919_pock.pdb --pocket_dual ../conf_B_4240_pock.pdb --score -7 --score_dual 0 --start 1 --qed 0.3 --processor 30
```

An example, for the paper data reproduce, single conformation A:
```
cd dual-conf-3D-MCTS/SARS-CoV-2/Mpro_A_B_conformation_generate_analysis/run_conf_A_dual_QWEN_test/
python ../../../3D-MCTS.py --num_sims 100000 --ligand ../ligand_A_s.sdf --protein ../conf_A_919.pdb --pocket ../conf_A_919_pock.pdb --score -7 --start 1 --frag_lib '../../../frags/fragments.txt' --qed 0.3 --processor 30
```

An example, for the paper data reproduce, single conformation B:
```
cd dual-conf-3D-MCTS/SARS-CoV-2/Mpro_A_B_conformation_generate_analysis/run_conf_B_dual_QWEN_test/
python ../../../3D-MCTS.py --num_sims 100000 --ligand ../ligand_B_s.sdf --protein ../conf_B_4240.pdb --pocket ../conf_B_4240_pock.pdb --score -7 --start 1 --frag_lib '../../../frags/fragments.txt' --qed 0.3 --processor 30
```



Reference
--------
1. Reference the single-target 3D-MCTS lib from: https: https://github.com/Brian-hongyan/3D-MCTS
2. Reference the gnina from: https://github.com/gnina/gnina
3. Reference the MPerformer lib from: https://github.com/FanmengWang/MPerformer
4. Reference the Uni-Core from: https://github.com/dptech-corp/Uni-Core/releases
5. Reference the ADFR from: https://ccsb.scripps.edu/adfr/


bash ./data/get_data.sh
python ./data/scripts/TransientDMD_surrogates-PyDMD_JOV_Latest.py
pdflatex -shell-escape ./paper/dmd_lra_surrogate_ans.tex
pdflatex -shell-escape ./paper/presentation/LRA_ANS_presentation.tex

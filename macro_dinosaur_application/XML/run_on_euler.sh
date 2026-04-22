#!/bin/bash

# Ensure directories exist
mkdir -p logs
mkdir -p results

for suffix in A B C D E
do

    # 1. No GLM
    sbatch --time=120:00:00 \
           --job-name=noGLM${suffix} \
           --output=logs/noGLM_${suffix}_%j.log \
           --wrap="echo 'Job start: ' \$(date);
                   START=\$(date +%s);
                   beast -D "suffix=${suffix}" -statefile results/Dinosaurs_noGLM_${suffix}.state -overwrite XML/Dinosaurs_noGLM.xml;
                   END=\$(date +%s);
                   echo 'Job end: ' \$(date);
                   ELAPSED=\$((END-START));
                   echo 'Elapsed time (seconds):' \$ELAPSED"

    # 2. GLM 3 p indicators
    sbatch --time=120:00:00 \
       --job-name=GLM3p${suffix} \
       --output=logs/GLM_3p_${suffix}_%j.log \
       --wrap="echo 'Job start: ' \$(date);
               START=\$(date +%s);
               beast -D "suffix=${suffix},n_pred=3,preds=collectionCounts,formationCounts,gridCellCounts" \
                -statefile results/Dinosaurs_GLM_3p_${suffix}.state -resume XML/Dinosaurs_GLM_sumprior.xml;
               END=\$(date +%s);
               echo 'Job end: ' \$(date);
               ELAPSED=\$((END-START));
               echo 'Elapsed time (seconds):' \$ELAPSED"


    # 3. GLM 3p indicators + error distribution class
    sbatch --time=120:00:00 \
       --job-name=GLM3perrordistr${suffix} \
       --output=logs/GLM_3p_errordistr_${suffix}_%j.log \
       --wrap="echo 'Job start: ' \$(date);
               START=\$(date +%s);
               $JAVA -jar GLMPrior.jar -D "suffix=${suffix},n_pred=3,preds=collectionCounts,formationCounts,gridCellCounts" -statefile results/Dinosaurs_GLM_3p_error_distr_${suffix}.state -overwrite XML/Dinosaurs_GLM_errordistr.xml;
               END=\$(date +%s);
               echo 'Job end: ' \$(date);
               ELAPSED=\$((END-START));
               echo 'Elapsed time (seconds):' \$ELAPSED"

    # 4. GLM 12 p
    sbatch --time=120:00:00 \
       --job-name=GLM12p${suffix} \
       --output=logs/GLM_12p_${suffix}_%j.log \
       --wrap="echo 'Job start: ' \$(date);
               START=\$(date +%s);
               beast -D "suffix=${suffix},n_pred=12,preds=collectionCounts,collectionCounts_r1,collectionCounts_r2,collectionCounts_rev,formationCounts,formationCounts_r1,formationCounts_r2,formationCounts_rev,gridCellCounts,gridCellCounts_r1,gridCellCounts_r2,gridCellCounts_rev" \
                -statefile results/Dinosaurs_GLM_12p_${suffix}.state -resume XML/Dinosaurs_GLM_sumprior.xml;
               END=\$(date +%s);
               echo 'Job end: ' \$(date);
               ELAPSED=\$((END-START));
               echo 'Elapsed time (seconds):' \$ELAPSED"

    
done




# # 12p indicators
#     sbatch --time=120:00:00 \
#        --job-name=GLM12perrordistr${suffix} \
#        --output=logs/GLM_12p_errordistr_${suffix}_%j.log \
#        --wrap="echo 'Job start: ' \$(date);
#                START=\$(date +%s);
#                $JAVA -jar GLMPrior.jar -D "suffix=${suffix},n_pred=12,preds=collectionCounts,collectionCounts_r1,collectionCounts_r2,collectionCounts_rev,formationCounts,formationCounts_r1,formationCounts_r2,formationCounts_rev,gridCellCounts,gridCellCounts_r1,gridCellCounts_r2,gridCellCounts_rev" -statefile results/Dinosaurs_GLM_12p_error_distr_${suffix}.state -overwrite XML/Dinosaurs_GLM_errordistr.xml;
#                END=\$(date +%s);
#                echo 'Job end: ' \$(date);
#                ELAPSED=\$((END-START));
#                echo 'Elapsed time (seconds):' \$ELAPSED"



    # # 2. GLMon 3p
    # sbatch --time=120:00:00 \
    #    --job-name=GLMon3p${suffix} \
    #    --output=logs/GLMon_3p_${suffix}_%j.log \
    #    --wrap="echo 'Job start: ' \$(date);
    #            START=\$(date +%s);
    #            beast -D "suffix=${suffix},n_pred=3,preds=collectionCounts,formationCounts,gridCellCounts" -statefile results/Dinosaurs_GLMon_3p_${suffix}.state -overwrite XML/Dinosaurs_GLMon.xml;
    #            END=\$(date +%s);
    #            echo 'Job end: ' \$(date);
    #            ELAPSED=\$((END-START));
    #            echo 'Elapsed time (seconds):' \$ELAPSED"

    # # 3. GLMon 3p + error
    # sbatch --time=120:00:00 \
    #    --job-name=GLMon3perror${suffix} \
    #    --output=logs/GLMon_3p_error_${suffix}_%j.log \
    #    --wrap="echo 'Job start: ' \$(date);
    #            START=\$(date +%s);
    #            beast -D "suffix=${suffix},n_pred=3,preds=collectionCounts,formationCounts,gridCellCounts" -statefile results/Dinosaurs_GLMon_3p_error_${suffix}.state -overwrite XML/Dinosaurs_GLMon_error.xml;
    #            END=\$(date +%s);
    #            echo 'Job end: ' \$(date);
    #            ELAPSED=\$((END-START));
    #            echo 'Elapsed time (seconds):' \$ELAPSED"

    # # 3b. GLMon 3p + error distribution class
    # sbatch --time=120:00:00 \
    #    --job-name=GLMon3plserrordistr${suffix} \
    #    --output=logs/GLMon_3pls_errordistr_${suffix}_%j.log \
    #    --wrap="echo 'Job start: ' \$(date);
    #            START=\$(date +%s);
    #            $JAVA -jar GLMPrior.jar -D "suffix=${suffix},n_pred=3ls,preds=collectionCounts,formationCounts,gridCellCounts" -statefile results/Dinosaurs_GLMon_3pls_errordistr_${suffix}.state -overwrite XML/Dinosaurs_GLMon_errordistr.xml;
    #            END=\$(date +%s);
    #            echo 'Job end: ' \$(date);
    #            ELAPSED=\$((END-START));
    #            echo 'Elapsed time (seconds):' \$ELAPSED"


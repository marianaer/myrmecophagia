#!/bin/bash
TRINOTATE_PATH=/home/escobar/bin/Trinotate
TRINITY_PATH=/home/escobar/bin/trinityrnaseq-v2.14.0
trinotate_report=$1
count_matrix=$2

# Generate map of feature id to annotated feture id. from the Trinotate report file
${TRINOTATE_PATH}/util/Trinotate_get_feature_name_encoding_attributes.pl ${trinotate_report} > annotation_feature_map.txt

# Integrate expression matrixes with functional annotations
${TRINITY_PATH}/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl ${count_matrix} annotation_feature_map.txt > ${count_matrix}_annotated

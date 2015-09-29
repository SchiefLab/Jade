#Jared Adolf-Bryfogle (Dunbrack Lab, Schief Lab)
#Script for coloring and selection for AHo renumbered antibodies. 
#Use combinations for fast selectiong and analysis:
#  Ex: show sticks, epitope
#      show sticks, CDRs
#      show lines, L1 L2 L3 
#      hide lines, CDRs


select antibody, chain L or chain H
select antigen, not (chain L or chain H)

#L1 Selection coloring
select L1, chain L and resid 24-42
select L1_epitope, ((L1) around 5) and not (antibody)
select L1_frame, ((L1) around 5) and not (antigen)
show cartoon, L1
color yellow, L1 and symbol c

#L2 Selection coloring
select L2, chain L and resid 57-72
select L2_epitope, (L2 expand 5) and not (antibody)
select L2_frame, (L2 expand 5) and not (antigen)
color orange, L2 and symbol c
show cartoon, L2

#L3 Selection coloring
select L3, chain L and resid 107-138
select L3_epitope, (L3 around 5) and not antibody
select L3_frame, (L3 around 5) and not antigen
color salmon, L3 and symbol c
show cartoon, L3

#H1 Selection coloring
select H1, chain H and resid 24-42
select H1_epitope, (H1 around 5) and not antibody
select H1_frame, (H1 around 5) and not antigen
color cyan, H1 and symbol c
show cartoon, H1

#H2 Selection coloring
select H2, chain H and resid 57-69
select H2_epitope, (H2 around 5) and not antibody
select H2_frame, (H2 around 5) and not antigen
color slate, H2 and symbol c
show cartoon, H2

#H3 Selection coloring
select H3, chain H and resid 107-138
select H3_epitope, (H3 around 5) and not antibody
select H3_frame, (H3 around 5) and not antigen
color magenta, H3 and symbol c
show cartoon, H3

#Groups
#group CDRs, L1 L2 L3 H1 H2 H3
#group epitopes, L1_epitope L2_epitope L3_epitope H1_epitope H2_epitope H3_epitope
#group frames, L1_frame L2_frame L3_frame H1_frame H2_frame H3_frame
select CDRs, L1 or L2 or L3 or H1 or H2 or H3
select epitope, L1_epitope or L2_epitope or L3_epitope or H1_epitope or H2_epitope or H3_epitope
select full_epitope, antibody around 5
select full_paratope, antigen around 5
group L1_group, L1 L1_epitope L1_frame
group L2_group, L2 L2_epitope L2_frame
group L3_group, L3 L3_epitope L3_frame
group H1_group, H1 H1_epitope H1_frame
group H2_group, H2 H2_epitope H2_frame
group H3_group, H3 H3_epitope H3_frame

#group full, full_epitope full_paratope






hide lines
hide (hydro)
util.cnc


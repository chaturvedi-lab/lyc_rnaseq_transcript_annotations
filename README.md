# Lycaeides RNA sequencing project UNR

#### Pre-requisiite

```
conda create -n dev python=3.8
conda activate dev

pip install pandas goatools tqdm
```

#### SNP Annotation
```
python create_snp_annotations.py --pos mappos_sub.txt --ann annot1631_sub.txt --out sub_snp_annotations_table_1631.out
```
- log `SnpPos(2600), ScaffId(1648): 100%|██████████████████████████| 20/20 [00:02<00:00,  9.52it/s]`

#### Transcript Annotation
```
python create_transcript_annotations.py --pos result1_all.txt --ann annot1631_sub.txt --out sub_transcript_annotations_table_1631.out
```
- log `CuffId(9069.1), ScaffId(1632): 100%|███████████████████████| 227/227 [00:24<00:00,  9.11it/s]`

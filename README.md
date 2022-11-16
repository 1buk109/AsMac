# AsMac

AsMac is a machine learning framework designed for pairwise comparison of [ribosomal RNA](https://en.wikipedia.org/wiki/Ribosomal_RNA) sequences.

## 1. Installation

- install python

```
pyenv install 3.10.6
```

- built new virtualenv

```
pyenv virtualenv 3.10.6 AsMac
```

- switch virtualenv

```
pyenv local 3.10.6/envs/AsMac
```

Build the cython executable by:

```
python setup_softnw.py build_ext --inplace
```

## 2. Usage

Use -h or --help flags to get help text for the program.<br />

```
python AsMac.py -h
```

Choose a model with the -m or --model flags, and run the script like the example:

```
python AsMac.py -i path/to/input.fasta -o path/to/output.csv -m 16S-full
```

AsMac takes input sequences written in fasta format and outputs a pairwise distance matrix in a CSV table. For example:
| | seq_1 | seq_2 | seq_3 |
| ----------- | ----------- | ----------- | ----------- |
| seq_1 | 0 | 0.2 | 0.3 |
| seq_2 | 0.2 | 0 | 0.1 |
| seq_3 | 0.3 | 0.1 | 0 |

## 3. Train new model (optional)

## Preparing training and testing data

Download `C++` library [SeqAn](https://github.com/seqan/seqan) or use the attached version to compile the c++ code for NW alignment.

```
unzip seqan.zip
g++ -I . -std=c++1z -o CppAlign/align CppAlign/main.cpp CppAlign/read_fasta.cpp
```

Then generate alignment distance results for the input sequences.

```
./CppAlign/align ./data/training_seq.fa 0
```

The alignment process is very time-consuming on a personal computer. The demo code uses the finished result: training_dist_prepared.txt

## Run the demo code for training

```
jupyter notebook demo.ipynb
```

## Reference

1. [Jian Chen, Le Yang, Lu Li, Steve Goodison, and Yijun Sun. "Alignment-free Comparison of Metagenomics Sequences via Approximate String Matching." Bioinformatics Advances (2022).](https://academic.oup.com/bioinformaticsadvances/advance-article/doi/10.1093/bioadv/vbac077/6767843)<br />
2. [Koide, Satoshi, Keisuke Kawano, and Takuro Kutsuna. "Neural edit operations for biological sequences." Advances in Neural Information Processing Systems 31 (2018).](https://proceedings.neurips.cc/paper/2018/hash/d0921d442ee91b896ad95059d13df618-Abstract.html)
3. https://github.com/puar-playground/AsMac

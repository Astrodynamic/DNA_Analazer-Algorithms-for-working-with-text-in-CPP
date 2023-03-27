# Classical bioinformatics problems 

Implementation of the DNA Analyzer project.


## Contents

1. [Chapter I](#chapter-i) \
    1.1. [Introduction](#introduction)
2. [Chapter II](#chapter-ii) \
    2.1. [Information](#information)
3. [Chapter III](#chapter-iii) \
    3.1. [Part 1](#part-1-implementation-of-the-exact-dna-search-project)  
	3.2. [Part 2](#part-2-implementation-of-the-nw-sequence-alignment-project)  
	3.3. [Part 3](#part-3-matching-regular-expressions)  
	3.4. [Part 4](#part-4-k-similar-strings)  
	3.5. [Part 5](#part-5-bonus-minimum-window-substring)  

## Chapter I 

![Classical bioinformatics problems](misc/images/DNA_Analyzer.JPG)

Once again, Eve got into a strange new project, which happens often in her company. The Applied Bioinformatics Research department needed someone who could write complex algorithms. They were working with a partner company in Britain on a project that would allow for better prescribing of medications for patients, and that required a ton of algorithms. Since most of them were already formalized and did not require much knowledge of biology or immersion in the subject matter of the project, it was possible to take an outside person to help. \
All the stars aligned for Eve. Her experience in algorithms came in handy, as did her previous interaction with the guys from Advanced Solutions Inc., those same British partners, as part of the long-suffering tasks on labyrinths.

No one in the Applied Bioinformatics Research department, despite Eve's expectations, was wearing white coats or carrying flasks with dangerous viruses that could lock the population of the entire planet at home for several years. Instead, the same programmers were sitting there, searching through long, incomprehensible strings on their computer screens. They were even only two floors above Eve. She had to move here for a while to speed up the communication process. Eve wasn't too upset about it, figuring she could pump up her algorithmic skills on this project one more time. However, she strongly regretted not being able to take her favorite chair with her...

## Introduction

In this project you will be introduced to the classical bioinformatics problems of substring search and sequence alignment, as well as other string processing algorithms.
 
## Chapter II

## Information

An important part of modern bioinformatics is the analysis of molecular sequences. Sequences are strings over an arbitrary fixed alphabet, such as the alphabet `{ A, C, G, T}` (DNA sequence). Here are examples of DNA sequences that are parts of the HIV virus:

```
GGTCTCTCTGGTTAGACCAGATCTGAGC
CTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTG
AAGCCTCAAT
```

### Substring search

In molecular sequence analysis, it is constantly necessary to compare sequences with each other in one way or another. For example, sometimes it is necessary to determine the _exact match_ of two strings, i.e. a character-by-character match of strings with the same size. A more difficult but realistic case is the substring search task, that is, finding all occurrences of a short string (_needle_) of length `m` in a large string (_haystack_) of length `n`, `m <= n`. Besides bioinformatics, the substring search in different variations is common in various _full-text search_ tasks.

>*Example:* \
>In the string "Madam, I'm Adam" the substring `am' occurs twice, at positions 4 and 13.

>*Example:* \
>The relatively short HIV virus sequence of approximately 10,000 nucleotides (see the full sequence [here](https://www.ncbi.nlm.nih.gov/nuccore/AF033819)) contains two occurrences of the string `AAGCCTCTCAAT` consisting of ten nucleotides.

Naive substring search solutions are effective when strings reach thousands of characters in size; for strings of millions of characters or more, and for a massive search of a large number of queries, brute-force algorithms may be too slow. Effective substring search algorithms are not entirely trivial. One simple but relatively effective algorithm is the Rabin-Karp algorithm. This algorithm is based on the idea of _hashing_ substrings of length `m`, that is, calculating some function that matches each string with some number (hash). 

A match of the hash for the substring from "haystack" with the "needle" hash potentially indicates that the substring did occur. The algorithm surpasses naive solutions for input data of sufficient size if the hash of some substring can be computed by knowing the hash of the substring that comes before it.

### Sequence alignment

Molecular sequences evolve and change over time. As a result of mutations and other evolutionary changes, the sequences of closely related organisms become evolutionarily distant from each other over time and become less and less similar. To compare non-perfectly matched sequences, exact string-search algorithms are not enough. To match such sequences, an _alignment_ is used. Alignment of two sequences is their recording one under the other, in which additional  "gaps" are entered . These gaps are arranged in such a way that the letters of the sequences under each other match. For example, the alignment of two sequences `GGGCGACACTCCACCATAGA` and `GGCGACACCCACCATACAT` (small pieces taken from the beginning of two versions of the genomes of the hepatitis C virus), might look like this:

```
GGGCGACACTCCACCATAGA-
|| |||||| |||||||| |
GG-CGACAC-CCACCATACAT
```

Note that two sequences can correspond to a large number of different alignments (with different numbers of matching positions). From a biological point of view, alignment helps to identify those positions in the sequences that are most likely to be homologous, i. e., having a common evolutionary origin.

For learning purposes, we can assume that the goal of alignment is to maximize the number of matching characters in an alignment and minimize the number of non-matching ones. The classic solution to this problem is the Needleman-Wunsch algorithm, proposed in 1970 and based on the idea of dynamic programming. This algorithm relies on a "similarity" function that estimates the cost ("score") of matching or mismatching characters written under each other, and maximizes the overall "similarity" of the alignment.

>*Example:* \
>The score of matching is 1, the score of mismatching is -1, the score of a gap is -2. Then the alignment from the previous example reaches a total score of 10: 17 matching characters, 1 mismatching, 3 gaps.

Advanced versions of this algorithm use similarity matrices, where the score of matches and mismatches can be different for different pairs of characters. 

## Chapter III

**General** instructions for all parts:

- The program must be developed in C++ language of C++17 standard 
- The program code must be located in the src folder
- When writing code it is necessary to follow Google Style
- Do not use outdated language constructs and libraries
- Provide a Makefile for building the program and tests (with targets all, clean, tests, app)
- Prepare full coverage of all functions/methods used in the implementation of each task with unit-tests
 
 
- The program must have a console interface

## Part 1. Implementation of the Exact DNA search project

Develop a program for a full-text search using the **Rabin-Karp algorithm**. \
The program takes *two* files as input. They contain sequences `a` and `b` of length `n <= 10000` and `m <= 100` respectively, `m <= n`. The output of the program is a list of positions of the string `a` at which `b` occurs in `a`. 
Input example: \
File `datasets/HIV-1_AF033819.3.txt` and a file with the following contents:
```
AAGCCTCAATAAAGCTT
```

Output example:
```
65 9150 9182
```

#### Execution time and memory consumption check

The Unix-like operating systems have a utility `/usr/bin/time` (not to be confused with the `time` command in `bash`). You can check the execution time and consumption of a program with the following command:
```
/usr/bin/time -v PROGRAM
```
where `PROGRAM` corresponds to the name of the executable file. 
The output of `/usr/bin/time` might look like this:
```
    Command being timed: "./nw"
    User time (seconds): 0.00
    System time (seconds): 0.00
    Percent of CPU this job got: 75%
    Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.00
    Average shared text size (kbytes): 0
    Average unshared data size (kbytes): 0
    Average stack size (kbytes): 0
    Average total size (kbytes): 0
    Maximum resident set size (kbytes): 3284
    Average resident set size (kbytes): 0
    Major (requiring I/O) page faults: 0
    Minor (reclaiming a frame) page faults: 137
    Voluntary context switches: 1
    Involuntary context switches: 0
    Swaps: 0
    File system inputs: 0
    File system outputs: 0
    Socket messages sent: 0
    Socket messages received: 0
    Signals delivered: 0
    Page size (bytes): 4096
    Exit status: 0
```
You are interested in time (`Elapsed (wall clock) time`) and memory (`Maximum resident set size (kbytes)`).

Maximum execution time: 1 sec \
Maximum memory consumption: 128 MB

## Part 2. Implementation of the NW sequence alignment project

### Part 2.1 Calculation of the optimal score

Develop a program to align two sequences over the alphabet `{A, C, G, T}`. \
The input of the program is a file with *three* strings. The first string contains three numbers -- the score of the match, mismatch, and gap. The next two strings are the sequences for the alignment. The output of the program is *one* number -- the value of the score for the best alignment. 

Input example:
```
1 -1 -2
GGGCGACACTCCACCATAGA
GGCGACACCCACCATACAT
```

Output example:
```
10
```

### Part 2.2 Recovering optimal alignment

Add to the program a recovery of the optimal alignment, for which the maximum score is reached. The output of the program is the value of maximum score, and under it is a record of two strings one under the other with gaps. Matching characters at the same positions are marked with a vertical line. 
Output example:
```
10
GGGCGACACTCCACCATAGA-
|| |||||| |||||||| |
GG-CGACAC-CCACCATACAT
```

#### Execution time and memory consumption check

See instructions for the [first task](#part-1-implementation-of-the-exact-dna-search-project).

Maximum execution time: 1 sec \
Maximum memory consumption: 128 MB

## Part 3. Matching regular expressions

Develop a program to check if a sequence over the alphabet `{A, C, G, T}` matches a regular expression. \
The input of the program is a file with *two* strings. The first string contains the sequence for which a match will be checked. The second string contains a pattern containing characters from the alphabet and the following characters:

- `.` -- corresponds to any single character in the alphabet;
- `?` -- corresponds to any single character in the alphabet or the absence of a character;
- `+` --  corresponds to zero or more repetitions of the previous element;
- `*` -- corresponds to any sequence of characters from the alphabet or to the absence of characters.

The output of the program is *True*/*False* - whether the given sequence matches the pattern.

Input example:
```
GGCGACACCCACCATACAT
G?G*AC+A*A.
```
Output example:
```
True
```

**Note: Do not use ready-made libraries for regular expressions, such as _regex_ or _PCRE_**, in this part.

## Part 4. K-similar strings

Strings s1 and s2 are k-similar (for some non-negative integer *k*) if we can swap the positions of two letters in s1 exactly *k* times so that the resulting string equals s2.
Develop a program to check the k-similarity of two sequences over the alphabet `{A, C, G, T}`. \
The input of the program is a file with *two* strings. The output of the program is the smallest *k* for which s1 and s2 are k-similar. If the strings are not anagrams, output an error message.

Input example:
```
GGCGACACC
AGCCGCGAC
```

Output example:
```
3
```

## Part 5. Bonus. Minimum window substring

A substring is a continuous sequence of characters within a string.

Develop a program for minimum window substring for a sequence over the alphabet `{A, C, G, T}`. \
The input of the program is a file with *two* strings s and t. The window substring of string s is the substring containing all characters in string t (including duplicates).
The output of the program is a window substring of minimum length. If there is no window substring, return an empty string.

Input example:
```
GGCGACACCCACCATACAT
TGT
```

Output example:
```
GACACCCACCATACAT
```

*Explanation:*

The resulting substring must contain the characters `T`, `G`, `T`. The following substrings are suitable:
```
GACACCCACCATACAT
CGACACCCACCATACAT
GCGACACCCACCATACAT
GGCGACACCCACCATACAT
```

The result is a substring of _minimum_ length, so we choose the first one.


💡 [Tap here](https://forms.yandex.ru/u/635ab18369387220075b7edd/) **to leave your feedback on the project**. Pedago Team really tries to make your educational experience better.

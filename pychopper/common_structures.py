# -*- coding: utf-8 -*-

from collections import namedtuple

Hit = namedtuple('Hit', 'Ref RefStart RefEnd Query QueryStart QueryEnd Score')
Seq = namedtuple('Seq', 'Id Name Seq Qual Umi')
Segment = namedtuple('Segment', 'Left Start End Right Strand Len')

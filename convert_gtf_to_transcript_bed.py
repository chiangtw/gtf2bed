#! /usr/bin/env python

import argparse
import gzip
import re
from collections import namedtuple
from operator import itemgetter
from itertools import groupby


class Bed:
    BED_TITLE = (
        'chr',
        'start',
        'end',
        'name',
        'score',
        'strand',
        'thickStart',
        'thickEnd',
        'itemRGB',
        'blockCount',
        'blockSizes',
        'blockStarts'
    )
    _region = namedtuple('Region', ('start', 'end'))

    def __init__(self, regions, union=False,
                 name='.', score=0, thickStart=None, thickEnd=None):
        self.regions = regions
        self.union = union

        self.name = name
        self.score = score
        self._thickStart = thickStart
        self._thickEnd = thickEnd

        if not self.union:
            self.regions = self.reverse_regions_if_minus_strand(self.regions)
        else:
            self.regions = self.get_union_regions(self.regions)

        self._parse_regions(self.regions)

    @staticmethod
    def reverse_regions_if_minus_strand(regions):
        strand = regions[0][3]
        if strand == '-':
            regions = list(reversed(regions))

        return regions

    def _parse_regions(self, regions):
        self._parse_region(regions[0])
        for region in regions[1:]:
            self._add_block(region)

    def _parse_region(self, region):
        self.chrom = region[0]
        self.start = int(region[1]) - 1
        self.end = int(region[2])
        self.strand = region[3]

        self.block_count = 1
        self.block_sizes = [self.end - self.start]
        self.block_starts = [0]

    def _add_block(self, region):
        other = Bed([region])

        assert self.strand == other.strand

        self.block_count += other.block_count
        self.block_sizes += other.block_sizes

        all_starts = [
            start + self.start for start in self.block_starts
        ] + [other.start]

        if other.start < self.start:
            self.start = other.start

        self.block_starts = [start - self.start for start in all_starts]

        if other.end > self.end:
            self.end = other.end

    def get_data(self, all_fields=False):
        fields = [
            self.chrom,
            self.start,
            self.end,
            self.name,
            self.score,
            self.strand
        ]

        if all_fields:
            fields += [
                self.thickStart,
                self.thickEnd,
                0,
                self.block_count,
                self._list_to_str(self.block_sizes, sep=','),
                self._list_to_str(self.block_starts, sep=',')
            ]

        return fields

    def to_string(self, all_fields=False):
        data = self.get_data(all_fields=all_fields)
        bed_txt = self._list_to_str(data)
        return bed_txt

    @staticmethod
    def _list_to_str(list_, sep='\t', end=''):
        string = sep.join(map(str, list_)) + end
        return string

    @classmethod
    def get_union_regions(cls, regions):
        assert len(set(map(itemgetter(0), regions))) == 1, \
            "Not all regions in the same chromosome!"
        assert len(set(map(itemgetter(3), regions))) == 1, \
            "Not all regions at the same strand!"

        chr_ = regions[0][0]
        strand = regions[0][3]

        regions = sorted(
            map(
                lambda r: cls._region(r[1], r[2]),
                regions
            ),
            key=lambda r: r.start
        )

        union_regions = []
        r1 = regions[0]
        for r2 in regions[1:]:
            if r1.end < r2.start:
                union_regions.append(r1)
                r1 = r2
            else:
                if r1.end < r2.end:
                    r1 = cls._region(r1.start, r2.end)
        else:
            union_regions.append(r1)

        union_regions = tuple((chr_, start, end, strand)
                              for start, end in union_regions)
        return union_regions

    @property
    def thickStart(self):
        if self._thickStart:
            return self._thickStart
        else:
            return 0

    @property
    def thickEnd(self):
        if self._thickEnd:
            return self._thickEnd
        else:
            return 0


class GTF:
    _attr_patter = re.compile(r'([^ ;]+) \"?([^;"]+)\"?;')

    def __init__(self, raw_gtf_line):
        self._raw_gtf_line = raw_gtf_line

        self._init()
        self._parse(self._raw_gtf_line)

    def _init(self):
        self.chr_ = ''
        self.feature = ''
        self.start = 0
        self.end = 0
        self.strand = ''
        self.attrs = {}

        self.gene_id = ''
        self.transcript_id = ''

    @classmethod
    def _parse_attrs_string(cls, attrs_string):
        attrs = dict(re.findall(cls._attr_patter, attrs_string))
        return attrs

    def _parse(self, raw_gtf_line):
        if not raw_gtf_line.startswith('#'):
            data = raw_gtf_line.rstrip('\n').split('\t')

            self.chr_ = data[0]
            self.feature = data[2]
            self.start = int(data[3])
            self.end = int(data[4])
            self.strand = data[6]
            self.attrs = self._parse_attrs_string(data[8])

            self.gene_id = self.attrs.get('gene_id', '')
            self.transcript_id = self.attrs.get('transcript_id', '')


def read_file(file_):
    if file_.endswith('.gz'):
        return gzip.open(file_, 'rt')
    else:
        return open(file_)


def read_gtf_file(gtf_file):
    for line in read_file(gtf_file):
        if not line.startswith('#'):
            gtf = GTF(line)
            yield gtf


def get_transcripts_data(gtf_file):
    for gtf in read_gtf_file(gtf_file):
        if (gtf.transcript_id != '') and (gtf.feature != 'gene'):
            yield gtf


def convert_gtf_to_transcript_bed(gtf_file, show_CDS=False):
    groupby_tid = groupby(
        get_transcripts_data(gtf_file),
        key=lambda gtf: gtf.transcript_id
    )

    for tid, tgp in groupby_tid:
        tgp = list(tgp)

        exons = filter(lambda gtf: gtf.feature == 'exon', tgp)
        exon_regions = [
            (exon.chr_, exon.start, exon.end, exon.strand)
            for exon in exons
        ]

        if show_CDS:
            CDSs = filter(lambda gtf: gtf.feature == 'CDS', tgp)
            CDSs_start_end = [(CDS.start, CDS.end) for CDS in CDSs]

            if len(CDSs_start_end) > 0:
                # coding transcript
                coding_region_start = min(CDS[0] for CDS in CDSs_start_end) - 1
                coding_region_end = max(CDS[1] for CDS in CDSs_start_end)
            else:
                # non-coding transcript
                exon_regions_start = min(region[1] for region in exon_regions) - 1
                coding_region_start = exon_regions_start
                coding_region_end = exon_regions_start
                
            bed = Bed(
                exon_regions,
                name=tid,
                thickStart=coding_region_start,
                thickEnd=coding_region_end,
                union=True
            )

        else:
            bed = Bed(exon_regions, name=tid, union=True)

        bed_text = bed.to_string(all_fields=True)

        yield bed_text


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    parser.add_argument('-o', '--out_bed_file', type=argparse.FileType('w'))
    parser.add_argument('-CDS', '--show-CDS-region', action="store_true")

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    try:
        results = convert_gtf_to_transcript_bed(
            args.gtf_file,
            show_CDS=args.show_CDS_region
        )

        for line in results:
            print(line, file=args.out_bed_file)

    except BrokenPipeError:
        pass

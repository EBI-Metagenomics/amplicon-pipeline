#!/usr/bin/env python3

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(add_help="Convert the mapseq otutable to biome")
    parser.add_argument(
        "--otu-table",
        required=True,
        help=(
            "The OTU table produced for the taxonomies found in"
            " the reference databases that was used with MAPseq."
        ),
    )
    parser.add_argument(
        "--query",
        required=True,
        help="The output from the MAPseq that assigns a taxonomy to a sequence.",
    )
    parser.add_argument(
        "--out-file", required=True, help="The file storing the tsv file."
    )
    parser.add_argument(
        "--krona", required=False, help="Output file name for the Krona text"
    )
    parser.add_argument(
        "--label",
        required=False,
        default="Unspecified",
        help="Label to add to the top of the outfile OTU table.",
    )
    parser.add_argument("--no-tax-id-file", required=True, help="TODO")
    parser.add_argument("--taxid", action="store_true")

    args = parser.parse_args()

    tax_counter = {}

    with open(args.query, "r") as query_fh:
        tax = ""
        for line in query_fh.readlines():
            if line.startswith("#"):
                # Skip counts
                continue
            # Pull out the fields that we need
            line = line.strip()
            fields = line.split("\t")
            
            if len(fields) < 14:
                tax = "Unclassified"
            else:
                if not fields[13]:
                    tax = "Unclassified"
                else:
                    tax = fields[13]
                    while tax.endswith("__"):
                        fds = tax.split(";")
                        fds = fds[:-1]
                        tax = ";".join(fds)

            if tax in tax_counter:
                tax_counter[tax]["count"] += 1
            else:
                tax_counter[tax] = {"count": 1, "otu": 1, "taxid": ""}

    with open(args.otu_table, "r") as otu_table_fh:
        csv_reader = csv.reader(otu_table_fh, delimiter="\t")
        next(csv_reader)
        for row in csv_reader:
            # Simple two column table, OTU code and taxonomy string.
            otu = row[0]
            tax = row[1]
            tax_id = row[2] if len(row) >= 3 else None
            # Have we seen this tax string? Store the OTU
            if tax in tax_counter:
                tax_counter[tax]["otu"] = otu
                if tax_id:
                    tax_counter[tax]["taxid"] = tax_id

    for tax in tax_counter.keys():
        if tax_counter[tax]["otu"] is None:
            raise ValueError(f"Fatal, |{tax}| has not got an OTU code assigned")

    output_header = []
    notaxidfile_header = []

    if args.taxid:
        if args.label.startswith("UNITE"):
            output_header = [
                "# OTU ID",
                args.label,
                "taxonomy",
            ]
        else:
            output_header = [
                "# OTU ID",
                args.label,
                "taxonomy",
                "taxid",
            ]
        notaxidfile_header = [
            "# OTU ID",
            args.label,
            "taxonomy",
        ]
    else:
        output_header = [
            "# OTU ID",
            args.label,
            "taxonomy",
        ]
        notaxidfile_header = [
            "# OTU ID",
            args.label,
            "taxonomy",
        ]

    with open(args.out_file, "w") as output_fh:
        output_fh.write("# Constructed from biom file\n")
        output_csv = csv.writer(output_fh, delimiter="\t")
        output_csv.writerow(output_header)
        for tax in sorted(tax_counter.keys()):
            row = [
                tax_counter[tax]["otu"],
                f"{float(tax_counter[tax]['count']):.1f}",
                tax,
            ]
            if args.taxid:
                row.append(tax_counter[tax]["taxid"])
            output_csv.writerow(row)

    if args.taxid:
        with open(args.no_tax_id_file, "w") as notaxidfile_fh:
            notaxidfile_fh.write("# Constructed from biom file\n")
            csv_writer = csv.writer(notaxidfile_fh, delimiter="\t")
            csv_writer.writerow(notaxidfile_header)
            for tax in sorted(tax_counter.keys()):
                csv_writer.writerow(
                    [
                        tax_counter[tax]["otu"],
                        f"{float(tax_counter[tax]['count']):.1f}",
                        tax,
                    ]
                )

    if args.krona:
        with open(args.krona, "w") as k:
            for tax in sorted(tax_counter.keys()):
                tax_mod = tax.replace(r"\D_\d{1}__", "")
                tax_mod = "\t".join(tax_mod.split(";"))
                k.write(f"{tax_counter[tax]['count']}\t{tax_mod}\n")


if __name__ == "__main__":
    main()
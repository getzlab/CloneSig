import argparse
from CloneSig import CloneSigAnalyzer

def main():
    parser = argparse.ArgumentParser(description="Run CloneSig analysis for identifying genes with convergent mutations.")
    parser.add_argument("input_maf", help="Path to the MAF file")
    parser.add_argument("genes_list", help="Path to the file containing genes to test")
    parser.add_argument("--step", choices=['permutations', 'genes', 'all'], 
                      default='all', help='Which step to run')
    args = parser.parse_args()

    analyzer = CloneSigAnalyzer(args.input_maf, args.genes_list)
    analyzer.run_analysis(args.step)

if __name__ == "__main__":
    main()


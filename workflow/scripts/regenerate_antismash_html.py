#!/usr/bin/env python3
"""
Regenerate full antiSMASH HTML output from filtered GenBank files using antiSMASH's API.
"""

import sys
import shutil
import json
import logging
from pathlib import Path
from typing import List

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


def regenerate_antismash_html(gbk_dir: Path, output_dir: Path, bin_id: str, 
                               original_json: Path = None) -> int:
    """
    Regenerate antiSMASH HTML from GenBank files using antiSMASH's rendering engine.
    
    Args:
        gbk_dir: Directory containing filtered GenBank files
        output_dir: Output directory for HTML and other files
        bin_id: Bin identifier
        original_json: Optional path to original antiSMASH JSON file
        
    Returns:
        Number of records processed
    """
    gbk_dir = Path(gbk_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find GenBank files
    gbk_files = sorted(gbk_dir.glob("*.gbk"))
    
    if not gbk_files:
        logging.warning(f"No GenBank files found in {gbk_dir}")
        create_empty_html(output_dir, bin_id)
        return 0
    
    try:
        # Import antiSMASH modules
        from antismash.common.secmet import Record
        from antismash.outputs.html import generator
        from antismash.outputs.html.generator import generate_html_sections
        import antismash.config
        from antismash.config.args import build_parser
        
        logging.info(f"Loading {len(gbk_files)} GenBank files...")
        
        # Copy GenBank files to output directory first
        for gbk_file in gbk_files:
            shutil.copy(gbk_file, output_dir / gbk_file.name)
        
        # Load records from GenBank files
        records = []
        for gbk_file in gbk_files:
            logging.info(f"Processing {gbk_file.name}")
            try:
                record = Record.from_genbank(str(gbk_file), taxon="bacteria")[0]
                records.append(record)
            except Exception as e:
                logging.error(f"Could not load {gbk_file}: {e}")
                continue
        
        if not records:
            logging.warning("No valid records loaded")
            create_empty_html(output_dir, bin_id)
            return 0
        
        logging.info(f"Successfully loaded {len(records)} records")
        
        # Create a minimal antiSMASH options object
        # We need to fake enough options for the HTML generator to work
        parser = build_parser()
        
        # Use minimal arguments - just what's needed for HTML generation
        args = [
            "--output-dir", str(output_dir),
            "--output-basename", bin_id,
            "--html-title", f"BGCs in Bin {bin_id}",
        ]
        
        # Parse options
        options = parser.parse_args(args)
        antismash.config.update_config(options)
        config = antismash.config.get_config()
        
        # Create results structure
        # The HTML generator needs this structure
        results = create_results_structure(records, original_json)
        
        # Generate HTML using antiSMASH's generator
        logging.info("Generating HTML sections...")
        
        try:
            # This is the main function that generates all HTML output
            from antismash.outputs.html import write
            
            # Write HTML output
            write.write(records, results, config)
            
            logging.info(f"Successfully generated antiSMASH HTML in {output_dir}")
            
        except Exception as e:
            logging.error(f"Error generating full HTML: {e}")
            logging.info("Falling back to simplified HTML generation")
            generate_simplified_html(records, output_dir, bin_id)
        
        return len(records)
        
    except ImportError as e:
        logging.error(f"Could not import antiSMASH modules: {e}")
        logging.info("Falling back to simple HTML")
        create_simple_html_fallback(gbk_files, output_dir, bin_id)
        return len(gbk_files)
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        create_simple_html_fallback(gbk_files, output_dir, bin_id)
        return len(gbk_files)


def create_results_structure(records: List, original_json: Path = None) -> dict:
    """
    Create a results structure compatible with antiSMASH's HTML generator.
    
    This mimics the structure that antiSMASH creates during analysis.
    """
    from antismash.common.module_results import ModuleResults
    
    results = {
        'records': [],
        'timings': {},
        'input_file': 'filtered_bin'
    }
    
    # If we have the original JSON, try to extract relevant results
    if original_json and original_json.exists():
        try:
            with open(original_json) as f:
                original_results = json.load(f)
                
            # Extract relevant parts of results for our filtered contigs
            # This is simplified - full implementation would match records to results
            results['timings'] = original_results.get('timings', {})
            
        except Exception as e:
            logging.warning(f"Could not load original JSON: {e}")
    
    # Add minimal results for each record
    for record in records:
        record_results = {
            'record_id': record.id,
            'modules': {}
        }
        results['records'].append(record_results)
    
    return results


def generate_simplified_html(records: List, output_dir: Path, bin_id: str):
    """
    Generate a simplified HTML using antiSMASH's basic rendering capabilities.
    """
    from antismash.outputs.html.generator import generate_html_sections
    
    try:
        # Generate basic HTML structure
        html_content = generate_basic_html_structure(records, bin_id)
        
        # Write index.html
        index_path = output_dir / "index.html"
        with open(index_path, 'w') as f:
            f.write(html_content)
            
        logging.info(f"Generated simplified HTML at {index_path}")
        
    except Exception as e:
        logging.error(f"Error in simplified HTML generation: {e}")
        create_simple_html_fallback([output_dir / f"{r.id}.gbk" for r in records], 
                                    output_dir, bin_id)


def generate_basic_html_structure(records: List, bin_id: str) -> str:
    """
    Generate basic HTML structure for antiSMASH-like visualization.
    """
    from Bio import SeqIO
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>antiSMASH results - Bin {bin_id}</title>
    <style>
        body {{
            font-family: 'Helvetica Neue', Arial, sans-serif;
            margin: 0;
            padding: 0;
            background-color: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 100%);
            color: white;
            padding: 30px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            margin: 0;
            font-size: 2.5em;
            font-weight: 300;
        }}
        .header p {{
            margin: 10px 0 0 0;
            opacity: 0.9;
        }}
        .container {{
            max-width: 1400px;
            margin: 30px auto;
            padding: 0 20px;
        }}
        .summary {{
            background: white;
            padding: 25px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            margin-bottom: 30px;
        }}
        .summary h2 {{
            margin-top: 0;
            color: #2a5298;
            border-bottom: 2px solid #2a5298;
            padding-bottom: 10px;
        }}
        .stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }}
        .stat {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .stat-value {{
            font-size: 2.5em;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        .stat-label {{
            opacity: 0.9;
            font-size: 0.9em;
        }}
        .record {{
            background: white;
            margin-bottom: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        .record-header {{
            background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
            color: white;
            padding: 20px;
            cursor: pointer;
        }}
        .record-header:hover {{
            opacity: 0.95;
        }}
        .record-header h3 {{
            margin: 0;
            font-size: 1.3em;
        }}
        .record-content {{
            padding: 25px;
        }}
        .region {{
            background: #f8f9fa;
            border-left: 4px solid #667eea;
            padding: 15px;
            margin: 15px 0;
            border-radius: 4px;
        }}
        .region-header {{
            font-weight: 600;
            color: #2a5298;
            margin-bottom: 10px;
        }}
        .product-tag {{
            display: inline-block;
            background: #667eea;
            color: white;
            padding: 5px 12px;
            border-radius: 15px;
            font-size: 0.85em;
            margin: 5px 5px 5px 0;
        }}
        .download-btn {{
            display: inline-block;
            background: #28a745;
            color: white;
            padding: 10px 20px;
            border-radius: 6px;
            text-decoration: none;
            margin-top: 10px;
            font-weight: 500;
        }}
        .download-btn:hover {{
            background: #218838;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #dee2e6;
        }}
        th {{
            background-color: #f8f9fa;
            font-weight: 600;
            color: #495057;
        }}
        tr:hover {{
            background-color: #f8f9fa;
        }}
    </style>
    <script>
        function toggleRecord(id) {{
            const content = document.getElementById('content-' + id);
            if (content.style.display === 'none') {{
                content.style.display = 'block';
            }} else {{
                content.style.display = 'none';
            }}
        }}
    </script>
</head>
<body>
    <div class="header">
        <h1>🧬 antiSMASH Results</h1>
        <p>Biosynthetic Gene Clusters in Bin {bin_id}</p>
    </div>
    
    <div class="container">
        <div class="summary">
            <h2>Summary</h2>
            <div class="stats">
                <div class="stat">
                    <div class="stat-value">{len(records)}</div>
                    <div class="stat-label">Contigs with BGCs</div>
                </div>
                <div class="stat">
                    <div class="stat-value">{sum(len([f for f in r.get_regions()]) for r in records)}</div>
                    <div class="stat-label">Total Regions</div>
                </div>
                <div class="stat">
                    <div class="stat-value">{sum(len([f for f in r.get_cds_features()]) for r in records)}</div>
                    <div class="stat-label">Total CDS</div>
                </div>
            </div>
        </div>
"""
    
    # Add each record
    for i, record in enumerate(records):
        regions = record.get_regions()
        cds_features = record.get_cds_features()
        
        html += f"""
        <div class="record">
            <div class="record-header" onclick="toggleRecord({i})">
                <h3>📄 {record.id}</h3>
                <p style="margin: 5px 0 0 0; opacity: 0.9;">
                    {len(record.seq):,} bp | {len(regions)} region(s) | {len(cds_features)} genes
                </p>
            </div>
            <div id="content-{i}" class="record-content">
"""
        
        # Add regions
        if regions:
            for region in regions:
                products = region.products
                product_str = ", ".join(products) if products else "Unknown"
                
                html += f"""
                <div class="region">
                    <div class="region-header">Region {region.get_region_number()}</div>
                    <p><strong>Location:</strong> {region.location.start:,} - {region.location.end:,} 
                       ({len(region.location):,} bp)</p>
                    <p><strong>Products:</strong> 
"""
                for product in products:
                    html += f'<span class="product-tag">{product}</span>'
                
                html += f"""
                    </p>
                    <p><strong>Genes in region:</strong> {len([f for f in region.cds_children])} CDS</p>
                </div>
"""
        else:
            html += "<p><em>No regions detected in this record</em></p>"
        
        # Add GenBank download link
        gbk_filename = f"{record.id}.gbk"
        html += f"""
                <a href="{gbk_filename}" class="download-btn" download>
                    📥 Download GenBank File
                </a>
            </div>
        </div>
"""
    
    html += """
    </div>
</body>
</html>
"""
    
    return html


def create_empty_html(output_dir: Path, bin_id: str):
    """Create HTML for bins with no BGCs."""
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>No BGCs - Bin {bin_id}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            display: flex;
            align-items: center;
            justify-content: center;
            margin: 0;
        }}
        .container {{
            background: white;
            padding: 50px;
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.3);
            text-align: center;
        }}
        h1 {{ color: #333; margin-bottom: 20px; }}
        p {{ color: #666; font-size: 1.1em; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🔍 No BGCs Found</h1>
        <p>Bin <strong>{bin_id}</strong> does not contain any biosynthetic gene clusters.</p>
    </div>
</body>
</html>
"""
    with open(output_dir / "index.html", 'w') as f:
        f.write(html_content)


def create_simple_html_fallback(gbk_files: List[Path], output_dir: Path, bin_id: str):
    """Fallback HTML when antiSMASH API fails completely."""
    from Bio import SeqIO
    
    bgc_data = []
    for gbk_file in gbk_files:
        try:
            for record in SeqIO.parse(gbk_file, "genbank"):
                regions = [f for f in record.features if f.type == "region"]
                for region in regions:
                    products = region.qualifiers.get("product", ["Unknown"])
                    bgc_data.append({
                        "contig": record.id,
                        "file": gbk_file.name,
                        "length": len(record),
                        "products": ", ".join(products) if isinstance(products, list) else products
                    })
        except Exception as e:
            logging.error(f"Error parsing {gbk_file}: {e}")
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>BGCs in Bin {bin_id}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 15px; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 20px; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:hover {{ background-color: #f5f5f5; }}
        .download {{ background-color: #27ae60; color: white; padding: 8px 15px; border-radius: 4px; text-decoration: none; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Biosynthetic Gene Clusters - Bin {bin_id}</h1>
        <p><strong>Total BGC regions:</strong> {len(bgc_data)}</p>
        <table>
            <tr><th>Contig</th><th>Products</th><th>Length (bp)</th><th>Download</th></tr>
"""
    
    for bgc in bgc_data:
        html += f"""
            <tr>
                <td>{bgc['contig']}</td>
                <td><em>{bgc['products']}</em></td>
                <td>{bgc['length']:,}</td>
                <td><a href="{bgc['file']}" class="download">Download</a></td>
            </tr>
"""
    
    html += """
        </table>
    </div>
</body>
</html>
"""
    
    with open(output_dir / "index.html", 'w') as f:
        f.write(html)


# Snakemake entry point
if __name__ == "__main__":
    try:
        # Access Snakemake variables
        gbk_dir = Path(snakemake.input.filtered_gbk)
        output_dir = Path(snakemake.params.output_dir)
        bin_id = snakemake.params.bin_id
        
        # Try to get original JSON if available
        original_json = None
        if hasattr(snakemake.input, 'antismash_json'):
            original_json = Path(snakemake.input.antismash_json)
        
        n_records = regenerate_antismash_html(gbk_dir, output_dir, bin_id, original_json)
        
        logging.info(f"Successfully processed {n_records} records for bin {bin_id}")
        
    except NameError:
        # Not running as Snakemake script
        import argparse
        
        parser = argparse.ArgumentParser()
        parser.add_argument("--gbk-dir", required=True)
        parser.add_argument("--output-dir", required=True)
        parser.add_argument("--bin-id", required=True)
        parser.add_argument("--json", help="Original antiSMASH JSON file")
        
        args = parser.parse_args()
        
        regenerate_antismash_html(
            Path(args.gbk_dir),
            Path(args.output_dir),
            args.bin_id,
            Path(args.json) if args.json else None
        )
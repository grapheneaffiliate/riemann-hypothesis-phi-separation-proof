#!/usr/bin/env python3
"""Convert readme.txt to HTML format."""

import re

def escape_html(text):
    """Escape HTML special characters."""
    return text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def convert_math(text):
    """Convert mathematical symbols to HTML entities."""
    replacements = {
        'ϕ': '&phi;', 'π': '&pi;', '∞': '&infin;', '√': '&radic;',
        '∑': '&sum;', '∏': '&prod;', '∫': '&int;', '≤': '&le;',
        '≥': '&ge;', '≠': '&ne;', '≈': '&asymp;', '⇒': '&rArr;',
        '⇔': '&hArr;', '∈': '&isin;', '⊂': '&sub;', '∀': '&forall;',
        '∃': '&exist;', 'ℝ': '&reals;', 'ℂ': '&#8450;', 'ℤ': '&#8484;',
        'ℚ': '&#8474;', 'Δ': '&Delta;', 'Γ': '&Gamma;', 'Θ': '&Theta;',
        'Λ': '&Lambda;', 'Σ': '&Sigma;', 'ξ': '&xi;', 'ζ': '&zeta;',
        'ρ': '&rho;', 'σ': '&sigma;', 'γ': '&gamma;', 'δ': '&delta;',
        'ε': '&epsilon;', 'τ': '&tau;', 'μ': '&mu;', 'λ': '&lambda;',
        '✓': '&check;', '✗': '&#10007;', '∎': '&#8718;', '│': '|',
        '├': '|--', '└': '`--', '▼': '&darr;', '⋯': '&hellip;',
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text

def main():
    # Read the text file
    with open('readme.txt', 'r', encoding='utf-8') as f:
        content = f.read()

    # HTML header with CSS
    html = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>The φ-Separation Proof of the Riemann Hypothesis</title>
    <style>
        body { 
            font-family: 'Georgia', serif; 
            max-width: 900px; 
            margin: 0 auto; 
            padding: 20px; 
            line-height: 1.7; 
            background: #fefefe; 
            color: #333; 
        }
        h1 { 
            color: #1a1a2e; 
            border-bottom: 3px solid #d4af37; 
            padding-bottom: 10px; 
            text-align: center;
        }
        h2 { 
            color: #16213e; 
            margin-top: 40px; 
            border-left: 4px solid #d4af37; 
            padding-left: 15px; 
            background: #f8f9fa;
            padding: 10px 15px;
        }
        h3 { color: #0f3460; margin-top: 25px; }
        h4 { color: #1a1a2e; }
        .abstract { 
            background: #f8f9fa; 
            padding: 20px; 
            border-left: 5px solid #d4af37; 
            margin: 20px 0; 
        }
        .theorem { 
            background: #e8f4f8; 
            padding: 15px 20px; 
            border: 1px solid #b8d4e3; 
            margin: 15px 0; 
            border-radius: 5px; 
        }
        .proof { 
            background: #f0f7e6; 
            padding: 15px 20px; 
            border: 1px solid #c5deb8; 
            margin: 15px 0; 
            border-radius: 5px; 
        }
        .conclusion { 
            background: #fff3cd; 
            padding: 20px; 
            border: 2px solid #d4af37; 
            margin: 20px 0; 
            border-radius: 8px; 
            text-align: center; 
        }
        .conclusion h2 { background: none; border: none; text-align: center; }
        pre { 
            background: #2d2d2d; 
            color: #f8f8f2; 
            padding: 15px; 
            overflow-x: auto; 
            border-radius: 5px; 
        }
        code { 
            background: #e9ecef; 
            padding: 2px 6px; 
            border-radius: 3px; 
            font-family: 'Consolas', monospace; 
        }
        table { 
            border-collapse: collapse; 
            width: 100%; 
            margin: 15px 0; 
        }
        th, td { 
            border: 1px solid #ddd; 
            padding: 10px; 
            text-align: left; 
        }
        th { background: #16213e; color: white; }
        tr:nth-child(even) { background: #f8f9fa; }
        blockquote { 
            border-left: 4px solid #d4af37; 
            margin: 20px 0; 
            padding: 10px 20px; 
            background: #f8f9fa; 
            font-style: italic; 
        }
        .qed { float: right; font-weight: bold; }
        a { color: #0066cc; }
        .part-header { 
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%); 
            color: white; 
            padding: 15px 20px; 
            margin: 30px 0 20px 0;
            border-radius: 5px;
        }
        .part-header h2 { 
            color: white; 
            margin: 0; 
            border: none;
            background: none;
            padding: 0;
        }
        .metadata { 
            color: #666; 
            font-size: 0.95em; 
            margin: 5px 0; 
            text-align: center;
        }
        ul, ol { margin: 10px 0; padding-left: 30px; }
        li { margin: 5px 0; }
        .step { margin: 10px 0; padding-left: 20px; }
        .remark { 
            background: #fff8e1; 
            padding: 15px; 
            border-left: 4px solid #ffc107; 
            margin: 15px 0; 
        }
        hr { border: none; border-top: 2px solid #d4af37; margin: 30px 0; }
    </style>
</head>
<body>
'''

    lines = content.split('\n')
    i = 0
    in_abstract = False
    in_theorem = False
    in_proof = False
    in_list = False
    
    while i < len(lines):
        line = lines[i].rstrip()
        
        # Skip empty lines
        if not line.strip():
            if in_list:
                html += '</ul>\n'
                in_list = False
            html += '<br>\n'
            i += 1
            continue
        
        processed_line = convert_math(escape_html(line))
        
        # Main title
        if 'Proof of the Riemann Hypothesis' in line and line.startswith('The'):
            html += f'<h1>{processed_line}</h1>\n'
        
        # Metadata
        elif line.startswith('Author:') or line.startswith('AI Collaborators:') or line.startswith('Date:') or line.startswith('Framework:'):
            html += f'<p class="metadata"><strong>{processed_line}</strong></p>\n'
        
        # Part headers
        elif line.startswith('PART ') or line.startswith('APPENDIX '):
            if in_abstract:
                html += '</div>\n'
                in_abstract = False
            html += f'<div class="part-header"><h2>{processed_line}</h2></div>\n'
        
        # Section headers
        elif line == 'ABSTRACT':
            html += '<h2>ABSTRACT</h2>\n<div class="abstract">\n'
            in_abstract = True
        
        elif line == 'MAIN THEOREM':
            if in_abstract:
                html += '</div>\n'
                in_abstract = False
            html += '<h2>MAIN THEOREM</h2>\n<div class="theorem">\n'
            in_theorem = True
        
        elif line == 'CONCLUSION':
            html += '<div class="conclusion"><h2>CONCLUSION</h2>\n'
        
        elif line == 'THE RIEMANN HYPOTHESIS IS TRUE':
            html += f'<h3 style="color: #d4af37; font-size: 1.5em;">{processed_line}</h3>\n'
        
        # Subsection headers (numbered like 1.1, 8.2, etc.)
        elif re.match(r'^\d+\.\d+\s+', line):
            html += f'<h3>{processed_line}</h3>\n'
        
        # Theorem/Corollary/Lemma headers
        elif line.startswith('Theorem') or line.startswith('Corollary') or line.startswith('Lemma') or line.startswith('Proposition'):
            html += f'<div class="theorem"><h4>{processed_line}</h4>\n'
            in_theorem = True
        
        # Proof headers
        elif line.startswith('Proof'):
            html += f'<div class="proof"><h4>{processed_line}</h4>\n'
            in_proof = True
        
        # QED
        elif '∎' in line:
            html += f'<p>{processed_line}</p>\n'
            if in_proof:
                html += '</div>\n'
                in_proof = False
            elif in_theorem:
                html += '</div>\n'
                in_theorem = False
        
        # Steps
        elif line.startswith('Step '):
            html += f'<p class="step"><strong>{processed_line}</strong></p>\n'
        
        # Remark
        elif line.startswith('Remark'):
            html += f'<div class="remark"><strong>{processed_line}</strong>\n'
        
        # List items
        elif line.startswith('- ') or line.startswith('• '):
            if not in_list:
                html += '<ul>\n'
                in_list = True
            html += f'<li>{convert_math(escape_html(line[2:]))}</li>\n'
        
        # Blockquotes
        elif line.startswith('"') or (line.startswith('---') and '---' in line[3:]):
            html += f'<blockquote>{processed_line}</blockquote>\n'
        
        # Citation block
        elif line.startswith('@misc{'):
            html += '<pre><code>'
            while i < len(lines) and '}' not in lines[i]:
                html += escape_html(lines[i]) + '\n'
                i += 1
            if i < len(lines):
                html += escape_html(lines[i]) + '\n'
            html += '</code></pre>\n'
        
        # Repository/URL lines
        elif line.startswith('Repository:') or line.startswith('http'):
            url = line.replace('Repository: ', '').strip()
            if 'http' in url:
                html += f'<p><strong>Repository:</strong> <a href="{url}">{url}</a></p>\n'
            else:
                html += f'<p>{processed_line}</p>\n'
        
        # Keywords/Subject lines
        elif line.startswith('Keywords:') or line.startswith('Subject:') or line.startswith('Related Work:'):
            html += f'<p><em>{processed_line}</em></p>\n'
        
        # Regular paragraphs
        else:
            if in_list:
                html += '</ul>\n'
                in_list = False
            html += f'<p>{processed_line}</p>\n'
        
        i += 1
    
    # Close any open tags
    if in_list:
        html += '</ul>\n'
    if in_abstract:
        html += '</div>\n'
    if in_theorem:
        html += '</div>\n'
    if in_proof:
        html += '</div>\n'
    
    html += '''
</body>
</html>'''

    # Write HTML file
    with open('readme.html', 'w', encoding='utf-8') as f:
        f.write(html)
    
    print('Successfully created readme.html!')

if __name__ == '__main__':
    main()

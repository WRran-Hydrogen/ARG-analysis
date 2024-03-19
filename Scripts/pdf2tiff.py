import os
import argparse
from pdf2image import convert_from_path

def convert_pdf_to_tiff(pdf_path, output_path):
    images = convert_from_path(pdf_path)
    for i, image in enumerate(images):
        # 每个页面都将保存为提供的文件名
        image.save(output_path, 'TIFF')

def main():
    parser = argparse.ArgumentParser(description="Convert PDF to TIFF")
    parser.add_argument("pdf_path", help="Path to the PDF file")
    parser.add_argument("output_path", help="Full path for the output TIFF file")
    args = parser.parse_args()

    convert_pdf_to_tiff(args.pdf_path, args.output_path)

if __name__ == "__main__":
    main()


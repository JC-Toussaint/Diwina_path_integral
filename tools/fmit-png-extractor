#!/usr/bin/env python3
import os

def chunk_parser(fileName,verbose):
    import png
    img = png.Reader(fileName)
    comments = ""
    someMetaData = {
        "Pixel width": "unknown",
        "Pixel height": "unknown",
        "Pixel unit": "unknown",
        "Value min": "unknown",
        "Value max": "unknown",
        "Value unit": "unknown"
    }
    for balise in img.chunks():    
        chunk_type = balise[0].decode('ascii')
        if (verbose):
            data_length = len(balise[1])
            if (data_length < 40):
                print(f"{chunk_type}: {balise[1]}")
            else:
                print(f"{chunk_type}: ({data_length} bytes)")
        match chunk_type:
            case 'tEXt':
                idx = balise[1].find(b'\x00')
                key = balise[1][0:idx].decode('ascii')
                value = balise[1][idx+1:].decode('ascii')
                someMetaData[key] = value

            case 'iTXt':
                idx = balise[1].find(b'\x00')
                key = balise[1][0:idx].decode('utf-8')
                if balise[1][idx+1] != 0:
                    print(f'Cannot decode compressed field "{key}"')
                else:
                    idx = balise[1].rfind(b'\x00')  # find last NUL seprator
                    value = balise[1][idx+1:].decode('utf-8')
                    someMetaData[key] = value

            case 'IHDR':
                # http://www.libpng.org/pub/png/spec/1.2/PNG-Chunks.html
                width = int.from_bytes(balise[1][0:4],byteorder='big')
                height = int.from_bytes(balise[1][4:8],byteorder='big')
                bitDepth = int.from_bytes(balise[1][8:9],byteorder='big')
                colorType = int.from_bytes(balise[1][9:10],byteorder='big')
                compressionMethod = int.from_bytes(balise[1][10:11],byteorder='big')
                filterMethod = int.from_bytes(balise[1][11:12],byteorder='big')
                interlaceMethod = int.from_bytes(balise[1][12:13],byteorder='big')
                comments += " (" + str(width) + "x" + str(height) + " px, " + str(bitDepth) + " bit, Grayscale)"
            
            case 'pCAL':
                # https://w3c.github.io/png/extensions/Overview.html#C.pCAL
                idx = balise[1].find(b'\x00')
                calibrationName = balise[1][0:idx].decode('ascii')
                idx += 1
                originalZero = int.from_bytes(balise[1][idx:idx+4],byteorder='big')
                idx += 4
                originalMax = int.from_bytes(balise[1][idx:idx+4],byteorder='big')
                idx += 4
                equationType = int.from_bytes(balise[1][idx:idx+1],byteorder='big')
                idx += 1
                numberOfParameters = int.from_bytes(balise[1][idx:idx+1],byteorder='big')
                idx += 1
                idx2 = balise[1][idx:].find(b'\x00')
                someMetaData["Value unit"] = balise[1][idx:idx+idx2].decode('ascii')
                if (numberOfParameters == 2) and (equationType == 0) and (originalZero == 0) and (originalMax == 65535):
                    idx = balise[1].rfind(b'\x00')
                    prm1 = float(balise[1][idx+1:].decode('ascii'))
                    idx2 = balise[1][0:idx].rfind(b'\x00')
                    prm0 = float(balise[1][idx2+1:idx].decode('ascii'))
                    someMetaData["Value min"] = prm0
                    someMetaData["Value max"] = prm0 + prm1
                else:
                    print("cannot handle calibration function" + calibrationName)

            case 'sCAL':
                # https://w3c.github.io/png/extensions/Overview.html#C.sCAL
                sCAL_unit = int.from_bytes(balise[1][0:1],byteorder='big')
                match sCAL_unit:
                    case 1: #meter
                        str_unit = "m"
                    case 2: #radian
                        str_unit = "rad"
                    case _: #undefined
                        str_unit = "arbitrary units"
                
                someMetaData["Pixel unit"] = str_unit
                idx = balise[1].find(b'\x00')
                someMetaData["Pixel width"] = balise[1][1:idx].decode('ascii')
                idx += 1
                someMetaData["Pixel height"] = balise[1][idx:].decode('ascii')

    return [bitDepth,colorType,someMetaData,comments]

def get_params():
    import argparse
    description = 'convert .png file to tiff and save metada in a yaml file'
    epilogue = '''
    input png file is read, pCAL and sCAL and tEXt chunks are read and formatted in an output yaml file and image is converted to tiff   
    '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fileName',type=str,help='PNG input file')
    parser.add_argument('-o','--outputFileName',type=str,help='optional output file name, conversion will respect output file name extension if possible')
    parser.add_argument('-v','--verbose', help='verbose mode: print all png chunk names and some of their content', action="store_true")
    parser.add_argument('--version',action='version',version= __version__,help='show the version number')
    args = parser.parse_args()
    return args

def writeMetadata(fileName,metadata,comments):
    with open(fileName,'w') as file:
        file.write(comments)
        for key in metadata:
            value = metadata[key]
            if '\n' in value:
                value = "  " + value.rstrip().replace("\n","\n  ")  # add indentation
                if key == "Simulation parameters":
                    file.write(f"{key}:\n{value}\n")  # write it as a YAML sub-object
                else:
                    file.write(f"{key}: |\n{value}\n")  # write it as a YAML multiline string
            else:
                file.write(f"{key}: {value}\n")
        file.close()

def main():
    args = get_params()

    if not(os.path.exists(args.fileName)):
        message = "file " + args.fileName + " does not exist."
    else:
        comments = "# Metadata from " + args.fileName
        [bitDepth,colorType,metaData,moreComments] = chunk_parser(args.fileName,args.verbose)
        idx = args.fileName.find(".png")
        ymlFileName = args.fileName[0:idx] + ".yml"
        comments += moreComments +'\n'
        writeMetadata(ymlFileName,metaData,comments)
        message = "Metadata written to " + ymlFileName

        if (colorType==0):
            import numpy as np
            from PIL import Image
            # that piece of code is a workaround, since pillow does not seem to handle correctly 16bits png format we use numpy
            if (bitDepth==16):
                img = Image.fromarray(np.array(Image.open(args.fileName)).astype("uint16"))
            else:
                img = Image.open(args.fileName)
            if args.outputFileName == None:
                outputFileName = args.fileName[0:idx] + ".tif"
                img.save(outputFileName,'TIFF')
                message += ", image written to " + outputFileName + "."
            else:
                outputFileName = args.outputFileName
                try:
                    img.save(outputFileName)
                    message += ", image written to " + outputFileName + "."
                except:
                    print("Error: cannot do that conversion.")
        else:
            message += "\nError: cannot handle that png: colorType=" + str(colorType) + " not a Grayscale."

    print(message)

__version__ = '0.0.0'
if __name__ == "__main__":
    main()


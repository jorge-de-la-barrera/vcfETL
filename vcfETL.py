import csv
import gzip

from vcfETLExceptions import *


class vcfETL:
    def __init__(self,
                 filename: str,  # Name of the VCF file
                 normalized: bool = True,  # VCF file is normalized
                 sample_names_translation: dict = None,  # Dictionary with the translation of the sample names
                 ):
        self.params = {
            'normalized': normalized
        }
        self.filename = filename
        self.__open_vcf_file()  # Open VCF file and parse header

    def __open_vcf_file(self):
        if self.filename.endswith('.gz') or self.filename.endswith('.bgz'):
            self.file = gzip.open(self.filename, 'rt')
        else:
            self.file = open(self.filename, 'r')
        self.file_start = self.file.tell()
        self.__parse_header()

    def __parse_header(self):
        """
        Extrae la información del encabezado del fichero VCF y la almacena en un diccionario.
        Actualiza el puntero al fichero para que apunte a la línea inmediatamente posterior a la cabecera.
        Extrae los nombres de las muestras y los guarda como un atributo de la clase.
        Elimina '#' en el nombre de la primera columna.
        """

        def translate_sample_names(sample_names):
            """
            Traduce los nombres de las muestras a un formato más legible.
            """
            # TODO: Implementar esta función para que traduzca los nombres de las muestras
            return sample_names

        def extract_csq_header(line):
            """
            Extraxt CSQ header

            :param line: line from the VCF file
            :return: named-array
            """
            idx_start = line.find('Allele')
            return line[idx_start:-2].split('|')

        self.header = {}
        for line in self.file:
            if line.startswith('##'):
                key, value = line.strip()[2:].split('=', 1)
                if value.startswith('<ID=CSQ'):
                    # Extract CSQ header from vcf
                    self.csq_header = extract_csq_header(value)
                self.header.setdefault(key, []).append(value)
            elif line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                columns[0] = columns[0].replace('#', '')
                self.column_names = columns[:9]
                self.sample_names = translate_sample_names(columns[9:])
                break

    def reset(self):
        """
        Restablece el puntero del fichero al inicio del cuerpo del VCF.
        """
        self.file.seek(self.file_start)
        for line in self.file:
            if not line.startswith('#'):
                break

    def parse_line(self):
        """
        Generador que devuelve una variante del fichero VCF cada vez que es invocado.
        Almacena la información de cada muestra en un diccionario con clave el nombre de la muestra y como valor
        un diccionario cuya clave es cada una de las etiquetas leídas en el campo FORMAT (columna 9).
        """

        def unzip_single_info_item(item):
            ''' Devuelve un dict para un item del campo INFO '''

            def int_or_float_or_str(x):
                try:
                    return int(x)
                except ValueError:
                    try:
                        return float(x)
                    except ValueError:
                        return str(x)

            tmp = item.split('=', 1)
            if len(tmp) == 1:
                res = {tmp[0]: True}
            else:
                # Si es el campo CSQ, lo desglosamos
                if tmp[0] == 'CSQ':
                    csq = []
                    for allele in tmp[1].split(','):
                        item_annotation = dict(zip(self.csq_header, allele.split('|')))
                        if item_annotation.get('Feature_type') == 'Transcript':
                            csq.append(item_annotation)
                    res = {tmp[0]: csq}
                else:
                    res = {tmp[0]: [int_or_float_or_str(v) for v in tmp[1].split(',')]}
            return res

        field_handlers = {
            'CHROM': lambda x: x,
            'POS': lambda x: int(x),
            'ID': lambda x: x.split(';') if x != '.' else [],
            'REF': lambda x: x.upper(),
            'ALT': lambda x: x.split(',') if x != '.' else [],
            'QUAL': lambda x: float(x) if x != '.' else None,
            'FILTER': lambda x: x.split(';') if x != 'PASS' and x != '.' else ['PASS'],
            'INFO': lambda x: {k: v for item in x.split(';') for k, v in
                               unzip_single_info_item(item).items()} if x != '.' else {},
        }

        for line in self.file:
            if not line.startswith('#'):
                variant_data = line.strip().split('\t')
                variant_dict = {col: handler(data) for col, data, handler in zip(self.column_names, variant_data[:9],
                                                                                 field_handlers.values())}

                if (self.params.get('normalized') and len(variant_dict['ALT']) > 1) or not self.params.get(
                        'normalized'):
                    raise vcfETLExceptions.NonNormalizedVCFError('VCF file is not normalized. Please, normalize it before use vcfETL.')

                # Agrupamos la información de las muestras bajo la clave SAMPLES
                variant_dict['SAMPLES'] = {}
                format_keys = variant_data[8].split(':')
                for sample_name, sample_data in zip(self.sample_names, variant_data[9:]):
                    sample_values = sample_data.split(':')
                    variant_dict['SAMPLES'][sample_name] = {key: value for key, value in
                                                            zip(format_keys, sample_values)}

                var_id = '_'.join(
                    [variant_dict['CHROM'], str(variant_dict['POS']), variant_dict['REF'], variant_dict['ALT'][0]])
                yield var_id, variant_dict

    def parse_file(self):
        """
        Devuelve un diccionario con todas las variantes del fichero VCF.
        """
        self.reset()
        return {var_id: variant for var_id, variant in self.parse_line()}

    # Methods to use VCFParser as a context manager
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

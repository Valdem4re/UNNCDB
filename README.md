# NNCDB

Прототип базы данных кристаллических структур для PostgreSQL, формируемой из cif-файлов.

## Использование

Создайте базу данных на сервере:
createdb <db_name> && psql <db_name> -f unncdb_schema.sql

## Сборка утилиты для заполнения базы данных
В директории parser выполните следующие команды:

- mkdir -p build
- cd build
- cmake ..
- cmake -- build .

## Использование утилиты для заполнения базы данных информацией из cif-файла:
cif_parser <input.cif> "\<conninfo\>"

Пример:

cif_parser ccdc-2023-dec.cif "host=localhost dbname=mydb user=myuser password=mypass"

#pragma once

#include <libpq-fe.h>
#include <string>
#include <stdexcept>
#include <memory>

using PGConnPointer = std::unique_ptr<PGconn, decltype(&PQfinish)>;
using PGResultPointer = std::unique_ptr<PGresult, decltype(&PQclear)>;

void PGConnect(const std::string& conninfo, PGConnPointer& conn);

void PGExec(PGConnPointer& conn, const std::string& sql);

void PGCopyBegin(PGConnPointer& conn, const std::string& sql);

void PGCopyEnd(PGConnPointer& conn);

void PGCopySend(PGConnPointer& conn, const std::string& line);
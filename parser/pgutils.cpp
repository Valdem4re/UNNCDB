#include "pgutils.hpp"

namespace {
using PGResultPointer = std::unique_ptr<PGresult, decltype(&PQclear)>;

PGResultPointer exec_checked(PGConnPointer& conn,
                             const std::string& sql,
                             ExecStatusType expected)
{
    PGresult* raw = PQexec(conn.get(), sql.c_str());
    if (!raw)
        throw std::runtime_error("PQexec returned null");

    PGResultPointer res(raw, PQclear);
    if (PQresultStatus(res.get()) != expected) {
        std::string err = PQerrorMessage(conn.get());
        throw std::runtime_error("SQL failed (" + sql + "): " + err);
    }
    return res;
}
} // namespace

void PGConnect(const std::string& conninfo, PGConnPointer& conn)
{
    conn.reset(PQconnectdb(conninfo.c_str()));
    if (!conn || PQstatus(conn.get()) != CONNECTION_OK) {
        std::string err = conn ? PQerrorMessage(conn.get()) : "null connection";
        conn.reset();
        throw std::runtime_error("Connection to database failed: " + err);
    }
}

void PGExec(PGConnPointer& conn, const std::string& sql)
{
    (void)exec_checked(conn, sql, PGRES_COMMAND_OK);
}

void PGCopyBegin(PGConnPointer& conn, const std::string& sql)
{
    (void)exec_checked(conn, sql, PGRES_COPY_IN);
}

void PGCopyEnd(PGConnPointer& conn)
{
    if (PQputCopyEnd(conn.get(), nullptr) != 1)
        throw std::runtime_error("PQputCopyEnd failed");

    PGresult* raw = nullptr;
    while ((raw = PQgetResult(conn.get())) != nullptr) {
        PGResultPointer res(raw, PQclear);
        if (PQresultStatus(res.get()) != PGRES_COMMAND_OK) {
            std::string err = PQerrorMessage(conn.get());
            throw std::runtime_error("COPY end failed: " + err);
        }
    }
}

void PGCopySend(PGConnPointer& conn, const std::string& line)
{
    int rc;
    do {
        rc = PQputCopyData(conn.get(), line.data(), (int)line.size());
        if (rc == -1) {
            std::string err = PQerrorMessage(conn.get());
            throw std::runtime_error("PQputCopyData failed: " + err);
        }
        if (rc == 0 && PQflush(conn.get()) == -1) {
            throw std::runtime_error("PQflush failed during COPY");
        }
    } while (rc == 0);
}